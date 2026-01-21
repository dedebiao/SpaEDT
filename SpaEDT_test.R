##test
devtools::install_github("dedebiao/SpaEDT")
library(SpaEDT)

# 设置随机种子确保可重复性
set.seed(42)

cell_data <- data.frame(
  cell_id = paste0("cell_", 1:2000),
  x = runif(2000, 0, 1000),  # X坐标
  y = runif(2000, 0, 1000),  # Y坐标
  type = sample(c("Tumor", "Immune", "Stroma"), 2000, 
                replace = TRUE, prob = c(0.4, 0.3, 0.3))
)

head(cell_data)

## 1. proportion
cell_data <- SpaEcoRatio(cell_data,
                         x = "x",
                         y = "y",
                         type = "type",
                         k=50)
head(cell_data)

## 2. 空间特征计算 ----
### 2.1 细胞密度特征 ----
cell_data <- SpaDensityMk(cell_data, radius = 20)
head(cell_data)

### 2.2 细胞类型空间关系 ----
# 计算到肿瘤细胞的接近度
target_cell <- SpaProxTarget(cell_data,type,"Tumor",50)
head(target_cell)

## 3. 肿瘤-免疫微环境分析 ----
### 基质阻断分析 ----

cell_data <- SpaBarrierQ(
  cell_data,
  type,
  tumor_type = "Tumor",
  immune_type = "Immune",
  stroma_type = "Stroma"
)
head(cell_data)

## 4. 空间图分析 ----
edges <- SpaGeoNet(cell_data, k = 6, max_dist = 50)

## 5. 空间社区特征 ----图计算
g <- SpaComLeiden(cell_data, edges, resolution = 0.5, n_iter = 10)

## 6. 可视化图聚类
plot_net <- SpaNetPlot(g, node_alpha = 0.7)
print(plot_net)