# 加载必要的包
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(scales)


setwd("H:\\毕业生\\贾诗雨\\manuscript\\cladistics\\R1\\analysis\\abba")




# 读取物种顺序文件和数据文件
species_order_file <- "all.namelist.txt"  # 包含物种顺序的文件
input_file <- "Hes.indABBA.avg.txt"  # 包含z_value和gamma的数据文件

# 读取物种顺序文件，假设每行包含一个物种名
species_order <- readLines(species_order_file)
species_order <- rev(species_order)
# 初始化一个对称矩阵，所有值初始化为0
species_pair_matrix_z <- matrix(0, nrow = length(species_order), ncol = length(species_order))
species_pair_matrix_d <- matrix(0, nrow = length(species_order), ncol = length(species_order))

# 设置行和列的名称
rownames(species_pair_matrix_z) <- species_order
colnames(species_pair_matrix_z) <- species_order
rownames(species_pair_matrix_d) <- species_order
colnames(species_pair_matrix_d) <- species_order

# 读取数据文件并更新矩阵
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(data)

# 更新矩阵中的值
for (i in 1:nrow(data)) {
  sp1 <- data$sp1[i]
  sp2 <- data$sp2[i]
  z_value <- data$z_value[i]
  d_value <- data$d_value[i]
  
  # 确保物种顺序符合species_order
  row_sp1 <- which(species_order == sp1)
  col_sp2 <- which(species_order == sp2)
  
  # 对于左下角的热图，更新z_value矩阵
  if (sp1 != sp2) {  # 排除同一物种的情况
    species_pair_matrix_z[row_sp1, col_sp2] <- z_value
    species_pair_matrix_z[col_sp2, row_sp1] <- z_value  # 对称矩阵
    
    # 对于右上角的热图，更新gamma矩阵
    species_pair_matrix_d[row_sp1, col_sp2] <- d_value
    species_pair_matrix_d[col_sp2, row_sp1] <- d_value  # 对称矩阵
  }
}

# 转换矩阵为长格式，用于绘图
melt_z <- melt(species_pair_matrix_z)
melt_d <- melt(species_pair_matrix_d)

head(melt_z)

# 绘制左下角热图 (z_value)

p1 <- ggplot(melt_z, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlBu")),
                       values = scales::rescale(c(0, 20, 100)), # 设置中间值为0
                       name = "z_value")  +
  labs(title = "z_value Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# 绘制右上角热图 (gamma)
p2 <- ggplot(melt_d, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = rev(brewer.pal(9, "RdYlBu")),
                       values = scales::rescale(c(0, 0.15, 0.25)), # 设置中间值为0
                       name = "D") +
  labs(title = "D_value Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2


p2 <-  ggplot(melt_d, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "grey88", size = 0.15) +            # 细网格线
  coord_fixed() +                                       # 方格等宽
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FBE5D6", "#F6B298", "#EB7B62", "#D64B3A", "#9C1515"),
    values  = rescale(c(0, 0.10, 0.20, 0.30, 0.40, 0.46)), # 根据你的数据上限微调
    limits  = c(0, max(melt_d$value, na.rm = TRUE)),
    oob     = squish,
    na.value = "grey80",                                # NA 灰色遮罩
    name = expression(italic(f)[d])                     # 颜色条标题
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),                        # 去掉坐标文字
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey90", fill = NA),
    legend.position = "right"
  )

p2


??scale_fill_viridis
  
ggsave("abba.z.pdf", p1, width = 10, height = 8, dpi = 300)
ggsave("abba.d.pdf", p2, width = 10, height = 8, dpi = 300)

