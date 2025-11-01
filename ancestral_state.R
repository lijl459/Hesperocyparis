library(phytools)
library(ape)
library(RRphylo)

set.seed(123)


#map growth rate
#median
getwd()
tree = read.tree("16sp_time_tree.tre")
plot(tree, cex = 0.8)
nodelabels()

tree <- rotate(tree, 19)
plot(tree, cex = 0.8)
nodelabels()

tree <- rotate(tree, 20)
plot(tree, cex = 0.8)
nodelabels()

tree <- rotate(tree, 21)
plot(tree, cex = 0.8)
nodelabels()

tree <- rotate(tree, 22)
plot(tree, cex = 0.8)
nodelabels()

tree <- rotate(tree, 23)
tree <- rotate(tree, 25)
tree <- rotate(tree, 26)

tree <- rotate(tree, 24)
plot(tree, cex = 0.8)
nodelabels()



plot(tree)



data = read.csv("clim_median.csv", header= T, row.names = 1)
head(data)
plot(tree)


library(RColorBrewer)

pal <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
pal_rev <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)


## 通用：按你的格式绘制 contMap 并保存 JPG
plot_contmap_jpg <- function(var_vec, var_name, main_title, legend_title,
                             filename,
                             pal = c("#2E5276","#A3D6CA","#F6EBB6","#CD6116","#A01D14")) {
  
  ## ---- 计算 Blomberg's K 与 P ----
  phylo_result <- phylosig(tree, var_vec, method = "K", test = TRUE)
  K_val <- round(phylo_result$K, 3)
  P_val <- signif(phylo_result$P, 3)
  
  ## ---- 构建连续性状映射 ----
  obj  <- contMap(tree, var_vec, plot = FALSE)
  cols <- colorRampPalette(pal)(5)
  obj2 <- setMap(obj, cols)
  
  ## ---- 打开输出设备 ----
  jpeg(filename, width = 2000, height = 2400, res = 300)
  on.exit(dev.off(), add = TRUE)
  
  par(oma = c(0, 0, 3, 0))
  par(mar = c(5, 4, 2, 2) + 0.1)
  plot(obj2, type = "phylogram", legend = FALSE, ann = FALSE)
  
  ## ---- 标题 ----
  mtext(main_title, side = 3, outer = TRUE, line = 1, cex = 1.2)
  
  ## ---- 添加图例 ----
  usr   <- par("usr")
  x_pos <- usr[1] + 0.02 * (usr[2] - usr[1])  # 2% 宽度处
  y_pos <- usr[3] + 0.03 * (usr[4] - usr[3])  # 3% 高度处
  par(xpd = TRUE)
  add.color.bar(
    leg   = 0.5 * max(nodeHeights(tree)),
    cols  = obj2$cols,
    title = legend_title,
    lims  = round(range(var_vec), 1),
    prompt = FALSE,
    x = x_pos, y = y_pos,
    width = 0.05
  )
  
  
  ## ---- 在右上角标注 K 和 P（斜体，两行） ----
  # 行距根据字符高度自动设置
  line_height <- strheight("M", cex = 1.3)
  
  # 第一行 K
  text(
    x = 0 * (usr[2] - usr[1]),
    y = 0.2 * (usr[4] - usr[3]),
    labels = bquote(italic(K) == .(K_val)),
    cex = 1.1,
    adj = 0
  )
  
  # 第二行 P（在第一行下移一个 line_height）
  text(
    x = 0 * (usr[2] - usr[1]),
    y = 0.16 * (usr[4] - usr[3]),
    labels = bquote(italic(P) == .(P_val)),
    cex = 1.1,
    adj = 0
  )
  
  ## ---- RRphylo: 计算速率 & 搜索shift ----
  source("RRphylo functions_v2.R")
  RR <- RRphylo(tree = tree, y = var_vec)
  SS <- search.shift(RR, status.type   = "clade",
                           auto.recognize = "yes",
                           covariate="FALSE",
                           test.single= "no")
  
  st <- SS$`shift test single clades`   # 命名向量：names 是节点编号
  # 想只圈显著：用下面这一行替换上一行
  # st <- SS$`shift test single clades`[SS$`shift test single clades` < 0.05]
  node_ids <- as.numeric(names(st))
  
  ## ---- 把这些节点用红色半透明实心圆圈出来 ----
  if (length(node_ids) > 0) {
    lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    xx <- lp$xx; yy <- lp$yy
    
    # 设置透明度参数（0 = 完全透明，1 = 不透明）
    alpha_fill <- 0.8   # 填充透明度
    alpha_border <- 0.8 # 边框透明度
    
    points(xx[node_ids], yy[node_ids],
           pch = 21,                             # 实心圆（带边框）
           cex = 3,                            # 圆的大小
           lwd = 1.8,                            # 边线粗细
           col = rgb(0.9, 0.25, 0.15, 0.9),
           bg  = rgb(0.95, 0.45, 0.25, 0.5))       # 填充颜色（红+透明）
 
  
  par(xpd = FALSE)
}
}




## ====== BIO2 ~ BIO19 ======
bio_1  <- as.matrix(data)[, 1]
bio_2  <- as.matrix(data)[, 2]
bio_3  <- as.matrix(data)[, 3]
bio_4  <- as.matrix(data)[, 4]
bio_5  <- as.matrix(data)[, 5]
bio_6  <- as.matrix(data)[, 6]
bio_7  <- as.matrix(data)[, 7]
bio_8  <- as.matrix(data)[, 8]
bio_9  <- as.matrix(data)[, 9]
bio_10 <- as.matrix(data)[,10]
bio_11 <- as.matrix(data)[,11]
bio_12 <- as.matrix(data)[,12]
bio_13 <- as.matrix(data)[,13]
bio_14 <- as.matrix(data)[,14]
bio_15 <- as.matrix(data)[,15]
bio_16 <- as.matrix(data)[,16]
bio_17 <- as.matrix(data)[,17]
bio_18 <- as.matrix(data)[,18]
bio_19 <- as.matrix(data)[,19]
ele <- as.matrix(data)[,20]
AI  <- as.matrix(data)[,21]


source("RRphylo functions_v2.R")
RR <- RRphylo(tree = tree, y = bio_19)
SS <- search.shift(RR,
                   status.type   = "clade",
                   auto.recognize = "yes",
                   covariate="FALSE",
                   test.single= "no")


## 建议的完整标题与图例单位（按常规 WorldClim 定义）
## 温度类（BIO2–BIO11）：你之前用的是 ×10 刻度，这里保持一致 “°C ×10”
## 降水类（BIO12–BIO19）：单位 “mm”
## BIO3（等温性）：百分比 “%”
## BIO4（温度季节性）：标准差×100（无量纲）写作 “SD ×100”

plot_contmap_jpg(bio_1, "BIO1",
                 "Ancestral reconstruction of annual mean temperature (BIO1)",
                 "BIO1 (°C ×10)", "BIO1_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_2,  "BIO2",
                 "Ancestral reconstruction of mean diurnal temperature range (BIO2)",
                 "BIO2 (°C ×10)", "BIO2_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_3,  "BIO3",
                 "Ancestral reconstruction of isothermality (BIO3)",
                 "BIO3 (%)", "BIO3_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_4,  "BIO4",
                 "Ancestral reconstruction of temperature seasonality (BIO4)",
                 "BIO4 (SD ×100)", "BIO4_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_5,  "BIO5",
                 "Ancestral reconstruction of max temperature of warmest month (BIO5)",
                 "BIO5 (°C ×10)", "BIO5_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_6,  "BIO6",
                 "Ancestral reconstruction of min temperature of coldest month (BIO6)",
                 "BIO6 (°C ×10)", "BIO6_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_7,  "BIO7",
                 "Ancestral reconstruction of temperature annual range (BIO7)",
                 "BIO7 (°C ×10)", "BIO7_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_8,  "BIO8",
                 "Ancestral reconstruction of mean temperature of wettest quarter (BIO8)",
                 "BIO8 (°C ×10)", "BIO8_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_9,  "BIO9",
                 "Ancestral reconstruction of mean temperature of driest quarter (BIO9)",
                 "BIO9 (°C ×10)", "BIO9_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_10, "BIO10",
                 "Ancestral reconstruction of mean temperature of warmest quarter (BIO10)",
                 "BIO10 (°C ×10)", "BIO10_contMap.jpg", pal = pal)

plot_contmap_jpg(bio_11, "BIO11",
                 "Ancestral reconstruction of mean temperature of coldest quarter (BIO11)",
                 "BIO11 (°C ×10)", "BIO11_contMap.jpg", pal = pal)



plot_contmap_jpg(bio_12, "BIO12",
                 "Ancestral reconstruction of annual precipitation (BIO12)",
                 "BIO12 (mm)", "BIO12_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_13, "BIO13",
                 "Ancestral reconstruction of precipitation of wettest month (BIO13)",
                 "BIO13 (mm)", "BIO13_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_14, "BIO14",
                 "Ancestral reconstruction of precipitation of driest month (BIO14)",
                 "BIO14 (mm)", "BIO14_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_15, "BIO15",
                 "Ancestral reconstruction of precipitation seasonality (BIO15)",
                 "BIO15 (CV)", "BIO15_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_16, "BIO16",
                 "Ancestral reconstruction of precipitation of wettest quarter (BIO16)",
                 "BIO16 (mm)", "BIO16_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_17, "BIO17",
                 "Ancestral reconstruction of precipitation of driest quarter (BIO17)",
                 "BIO17 (mm)", "BIO17_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_18, "BIO18",
                 "Ancestral reconstruction of precipitation of warmest quarter (BIO18)",
                 "BIO18 (mm)", "BIO18_contMap.jpg", pal = pal_rev)

plot_contmap_jpg(bio_19, "BIO19",
                 "Ancestral reconstruction of precipitation of coldest quarter (BIO19)",
                 "BIO19 (mm)", "BIO19_contMap.jpg",
                 pal = pal_rev)

## ====== ele 与 AI ======
ele <- as.matrix(data)[,20]
AI  <- as.matrix(data)[,21]

plot_contmap_jpg(ele, "ele",
                 "Ancestral reconstruction of elevation (ele)",
                 "Elevation (m)", "ele_contMap.jpg",
                 pal = pal)

plot_contmap_jpg(AI, "AI",
                 "Ancestral reconstruction of aridity index (AI)",
                 "AI (index)", "AI_contMap.jpg",
                 pal = pal_rev)



