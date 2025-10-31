# mcmc_compare.R
# 读取多个 MCMCtree 的 mcmc.txt，计算各节点 Mean 与 95% HPD，并画对比图
# 需要的包：tidyverse, readr, stringr, purrr, coda, ggrepel(可选)
# install.packages(c("tidyverse","coda","ggrepel"))


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(coda)
  library(ggplot2)
  library(ggrepel)
  library(phytools)
  library(ape)
})

library(phytools)
library(ape)
library(RRphylo)

tr <- read.tree("timetree.name.tre")
plot(tr, cex = 0.5)

# 翻转整个树的分支顺序
tr_rev <- ape::rotateConstr(tr, rev(tr$tip.label))

plot(tr_rev, cex = 0.5)
nodelabels()


# ======== 用户区：填写你的文件路径与标签 ========
inputs <- c(
  "1.mcmc.txt",
  "2.mcmc.txt",
  "3.mcmc.txt",
  "4.mcmc.txt",
  "5.mcmc.txt"
)

labels <- c( # 与 inputs 一一对应；若留空将自动用上级目录名
  "SN(1.34, 0.0134, 0)",
  "SN(1.34, 0.153, 0)",
  "B(0.78–1.526)",
  "B(1.04–1.64)",
  "Gamma"
)

outdir <- "results_mcmc"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ======== 工具函数 ========

guess_label <- function(path) {
  # 用上级目录名作为标签
  b <- basename(dirname(normalizePath(path, mustWork = FALSE)))
  if (b == "" || b == "." || b == "/") return(basename(path))
  b
}

read_mcmc <- function(path) {
  # 读取 mcmc.txt，跳过以 # 开头的注释，自动空白分隔
  suppressMessages(
    read_table(path, comment = "#", col_types = cols(.default = col_guess()))
  )
}

pick_time_cols <- function(df) {
  # 选择节点时间列：常见命名 t1 / t_1 / t(1) / t_n1 等
  cand <- names(df)
  patt <- "^(t\\d+|t_\\d+|t\\(\\d+\\)|t_n\\d+)$"
  cols <- cand[str_detect(cand, patt)]
  if (length(cols) == 0) {
    # 兜底：所有以 t 开头且为数值的列
    num_cols <- cand[sapply(df, is.numeric)]
    cols <- num_cols[str_starts(num_cols, "t")]
  }
  cols
}

normalize_node_name <- function(x) {
  # 将 t(12)、t_12、t12、t_n12 统一成 t12
  x <- str_replace_all(x, "^t_n", "t")
  x <- str_replace_all(x, "^t_", "t")
  x <- str_replace_all(x, "^t\\(", "t")
  x <- str_replace_all(x, "\\)$", "")
  x
}

hpd95 <- function(x) {
  # 计算 95% HPD；若 coda 报错则退化为分位数
  x <- x[is.finite(x)]
  if (length(x) < 10) return(c(NA, NA))
  m <- try(HPDinterval(as.mcmc(x), prob = 0.95), silent = TRUE)
  if (inherits(m, "try-error")) {
    q <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
    return(as.numeric(q))
  } else {
    return(as.numeric(m[1, c("lower", "upper")]))
  }
}

summarize_one_run <- function(path, label = NULL) {
  df <- read_mcmc(path)
  tcols <- pick_time_cols(df)
  if (length(tcols) == 0) stop("没有识别到节点时间列：", path)
  
  df_t <- df %>% select(all_of(tcols))
  # 逐列（节点）计算 Mean/HPD
  res <- map_dfr(names(df_t), function(cn) {
    x <- df_t[[cn]]
    mn <- mean(x, na.rm = TRUE)
    hp <- hpd95(x)
    tibble(NodeRaw = cn, Mean = mn, HPD_L = hp[1], HPD_U = hp[2])
  })
  
  res %>%
    mutate(Node = normalize_node_name(NodeRaw),
           Run  = ifelse(is.null(label) || label == "", guess_label(path), label))
}

# ======== 主流程：汇总五次运行 ========

if (length(labels) != length(inputs)) {
  labels <- rep(NA_character_, length(inputs))
}

sum_list <- map2(inputs, labels, summarize_one_run)
sum_all  <- bind_rows(sum_list) %>%
  distinct(Node, Run, .keep_all = TRUE) %>%
  arrange(as.numeric(str_extract(Node, "\\d+")), Run)

# 保存统计表
write_csv(sum_all, file.path(outdir, "node_times_summary.csv"))

# ======== 画图 1：节点×方法 的 点+误差线（横向） ========

# 节点排序（按编号自然排序）
node_levels <- sum_all %>% pull(Node) %>% unique() %>%
  { .[order(as.numeric(str_extract(., "\\d+")))] }

p1 <- ggplot(sum_all, aes(x = Mean, y = factor(Node, levels = node_levels), color = Run)) +
  geom_errorbarh(aes(xmin = HPD_L, xmax = HPD_U), height = 0, alpha = 0.6) +
  geom_point(size = 2) +
  labs(x = "Divergence time (Ma)", y = "Node", color = "Root prior",
       title = "Comparison of node ages under different root calibrations (MCMCtree)",
       subtitle = "Points = posterior mean; bars = 95% HPD") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(outdir, "times_point_error.pdf"), p1, width = 8.5, height = 10)
ggsave(file.path(outdir, "times_point_error.png"), p1, width = 8.5, height =10, dpi = 300)

# ======== 画图 2：相对偏移（各节点对自身5次的均值归一化） ========
sum_rel <- sum_all %>%
  group_by(Node) %>%
  mutate(Relative = Mean / mean(Mean, na.rm = TRUE)) %>%
  ungroup()

p2 <- ggplot(sum_rel, aes(x = Run, y = Relative, group = Node, color = Node)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Root prior", y = "Relative time (Mean / Node-mean)",
       title = "Relative stability of node ages across root priors") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(min(0.8, min(sum_rel$Relative, na.rm = TRUE)),
                           max(1.2, max(sum_rel$Relative, na.rm = TRUE)))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "times_relative_lines.pdf"), p2, width = 8.5, height = 5.2)
ggsave(file.path(outdir, "times_relative_lines.png"), p2, width = 8.5, height = 5.2, dpi = 300)

# ======== 画图 3：热图（可选，展示各节点在不同先验下的均值） ========
mat <- sum_all %>%
  select(Node, Run, Mean) %>%
  pivot_wider(names_from = Run, values_from = Mean)

mat_long <- mat %>%
  pivot_longer(-Node, names_to = "Run", values_to = "Mean")

p3 <- ggplot(mat_long, aes(x = Run, y = factor(Node, levels = node_levels), fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Root prior", y = "Node", fill = "Mean (Ma)",
       title = "Heatmap of node age means under different root priors") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "times_heatmap.pdf"), p3, width = 7.8, height = 6.5)
ggsave(file.path(outdir, "times_heatmap.png"), p3, width = 7.8, height = 6.5, dpi = 300)

message("Done. 输出在: ", normalizePath(outdir))



# ======== 新增：将“全部节点(左)” 与 “48–61(右)”并排输出 ========
# 需要 patchwork 来拼图
# install.packages("patchwork")
library(patchwork)

# 1) 准备 48–61 的数据与顺序
df_in <- sum_all %>%
  mutate(NodeNum = as.numeric(stringr::str_extract(Node, "\\d+"))) %>%
  filter(NodeNum >= 48 & NodeNum <= 61)

node_levels_in <- df_in %>% pull(Node) %>% unique() %>%
  { .[order(as.numeric(stringr::str_extract(., "\\d+")))] }

# 2) Point + Error（均值+95%HPD）——右侧子图
p1_in <- ggplot(df_in, aes(x = Mean, y = factor(Node, levels = node_levels_in), color = Run)) +
  geom_errorbarh(aes(xmin = HPD_L, xmax = HPD_U), height = 0, alpha = 0.6) +
  geom_point(size = 2) +
  labs(x = "Divergence time (Ma)", y = "Node (48–61)", color = "Root prior",
       title = "Nodes 48–61 (subset)",
       subtitle = "Points = posterior mean; bars = 95% HPD") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# 拼：左=全部节点 p1，右=48–61 p1_in
p1_bi <- p1 + p1_in + plot_layout(widths = c(1, 1))
ggsave(file.path(outdir, "times_point_error_all_vs_48_61.pdf"), p1_bi, width = 12, height = 6.5)
ggsave(file.path(outdir, "times_point_error_all_vs_48_61.png"), p1_bi, width = 12, height = 6.5, dpi = 300)

# 3) Relative line（相对偏移）——右侧子图
sum_rel_in <- sum_all %>%
  mutate(NodeNum = as.numeric(stringr::str_extract(Node, "\\d+"))) %>%
  filter(NodeNum >= 48 & NodeNum <= 61) %>%
  group_by(Node) %>%
  mutate(Relative = Mean / mean(Mean, na.rm = TRUE)) %>%
  ungroup()

p2_in <- ggplot(sum_rel_in, aes(x = Run, y = Relative, group = Node, color = Node)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Root prior", y = "Relative time (Mean / Node-mean)",
       title = "Relative stability (Nodes 48–61)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(
    min(0.9, min(sum_rel_in$Relative, na.rm = TRUE)),
    max(1.1, max(sum_rel_in$Relative, na.rm = TRUE))
  )) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 拼：左=全部节点 p2，右=48–61 p2_in
p2_bi <- p2 + p2_in + plot_layout(widths = c(1, 1))
ggsave(file.path(outdir, "times_relative_lines_all_vs_48_61.pdf"), p2_bi, width = 12, height = 5.2)
ggsave(file.path(outdir, "times_relative_lines_all_vs_48_61.png"), p2_bi, width = 12, height = 5.2, dpi = 300)

# 4) Heatmap（均值热图）——右侧子图
mat_long_in <- sum_all %>%
  mutate(NodeNum = as.numeric(stringr::str_extract(Node, "\\d+"))) %>%
  filter(NodeNum >= 48 & NodeNum <= 61) %>%
  select(Node, Run, Mean)

p3_in <- ggplot(mat_long_in, aes(x = Run, y = factor(Node, levels = node_levels_in), fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Root prior", y = "Node (48–61)", fill = "Mean (Ma)",
       title = "Heatmap (Nodes 48–61)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 拼：左=全部节点 p3，右=48–61 p3_in
p3_bi <- p3 + p3_in + plot_layout(widths = c(1, 1))
ggsave(file.path(outdir, "times_heatmap_all_vs_48_61.pdf"), p3_bi, width = 12, height = 6.5)
ggsave(file.path(outdir, "times_heatmap_all_vs_48_61.png"), p3_bi, width = 12, height = 6.5, dpi = 300)


