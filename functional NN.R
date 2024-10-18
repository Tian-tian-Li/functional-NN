
library(fda)
library(Rcpp)
library(RcppEigen)
library(FRegSigCom)
library(refund)
library(MASS)
library(readxl)
library(caret)
library(DALEX)
library(shapviz)
library(reshape2)
library(ggplot2)

source("new_fos.universal.nonlinear.R")


# read data --------------------------------------------------------------------

#death probability data
city_male_deathdata_original <- read_excel("31个省--城镇--男性--死亡概率.xlsx",range = "B2:AF22") 
city_female_deathdata_original <- read_excel("31个省--城镇--女性--死亡概率.xlsx", range = "B2:AF22") 
rural_male_deathdata_original <- read_excel("31个省--农村--男性--死亡概率.xlsx", range = "B2:AF22") 
rural_female_deathdata_original <- read_excel("31个省--农村--女性--死亡概率.xlsx", range = "B2:AF22") 


# 创建年龄组向量
age_groups <- c("0岁", "1-4岁", "5-9岁", "10-14岁", "15-19岁", "20-24岁", "25-29岁", "30-34岁", 
                "35-39岁", "40-44岁", "45-49岁", "50-54岁", "55-59岁", "60-64岁", "65-69岁", 
                "70-74岁", "75-79岁", "80-84岁", "85-89岁", "90-94岁")

# 创建一个数据框
age_groups_df <- data.frame(
  number = seq_along(age_groups), # 生成从1到年龄组长度的序列
  group = age_groups # 使用年龄组向量作为组列
)


Y <- cbind(city_male_deathdata_original,city_female_deathdata_original,rural_male_deathdata_original,rural_female_deathdata_original)
Y <- t(as.matrix(Y))
colnames(Y) <- age_groups
Y <- Y[,10:20]

dim(Y)[2]  # 检查列数

age_groups_Y <- data.frame(
  number = seq_along(colnames(Y)), # 生成从1到年龄组长度的序列
  group = colnames(Y) # 使用年龄组向量作为组列
)


#指标数据
X <- read_excel("newdata_2.xlsx",col_names=TRUE,range="A1:AQ125")
X <- X[,-1] #删掉第一列


# extracting rows 
vec <- c(
  "性别", "城乡","平均受教育年限",  "人均GDP",
  "第三产业占比", "每千人口卫生技术人员数",
  "每十万人口高等教育平均在校学生数", "交通事故发生数",
  "各地区城乡居民和职工基本医疗保险情况",  "医疗保健支出占消费性支出",
  "极端低温天数", "SO2", "NO2",   "PM2.5",
  "年平均相对湿度", "人均用水量", "地质灾害", "突发环境事件次数", "森林覆盖率",
  "城市污水处理率",  "年平均气温", "人均水资源量", "农作物受灾面积合计"
)

X <- X[ ,vec]


# 添加随机种子
# set.seed(123)

ntot =dim(X)[1]
ntrain=100


# train model -------------------------------------------------------------

#初始化参数，训练模型的次数
repeat_model <- 100
number_agegroup <- dim(Y)[2]

# 初始化一个空列表，用于存储每次循环中的模型结果
fit_cv_list <- vector("list", repeat_model)

# 初始化一个列表来保存每次循环的I结果
I_results <- list()

# 初始化一个列表来保存每次循环的训练集和测试集
train_test_sets <- list()


start_time <- proc.time()

for (n.repeat_model in 1:repeat_model) {
  
  I=sample(1:ntot,ntrain)
  I_results[[n.repeat_model]] <- I
  
  X.train <- X[I,]
  Y.train <- Y[I, ]
  
  X.test <- X[-I,]
  Y.test <- Y[-I,]
  
  # 将训练集和测试集保存到列表中
  train_test_sets[[n.repeat_model]] <- list(train = list(X = X.train, Y = Y.train),
                                            test = list(X = X.test, Y = Y.test))
  
  n.u <- NULL
  fit.cv <- cv.fos.universal.nonlinear(X.train, Y.train, t.n.basis = 6, u.n.basis = 6, n.u = n.u)
  fit_cv_list[[n.repeat_model]] <- fit.cv
  
}


end_time <- proc.time()
elapsed_time <- end_time - start_time
print(elapsed_time)



# 创建解释器----

explainer_non.fos_agegroup_list <- vector("list", number_agegroup)

# 根据 choose_age 的值来创建不同的 explainer_non.fos
for (i in 1:dim(Y)[2]) {
  
  # col_name <- colnames(Y)[i]
  # number_j <- age_groups_df$number[age_groups_df$group == col_name] #注意！！！
  
  for(n.repeat_model in 1:repeat_model){
    fit.cv <- fit_cv_list[[n.repeat_model]]
    explainer_non.fos_agegroup_list[[i]][[n.repeat_model]] <- explain(
      model = fit.cv,
      data = X,
      y = Y[, i],
      predict_function = get(paste0("pred.fos.universal.nonlinear.", i)),
      label = paste0("My Model.", colnames(Y)[i])
    ) 
  }
}




# 置换重要性---------------------------------------------------------------

# 初始化一个空列表，用于存储每次循环中的 explainer 结果
#repeat_perm <- 1 #在函数中设置即可。不用单独算了。


# * 分年龄段----
# 初始化一个空的列表，用于存储当前模型重复的年龄组实验结果

start_time <- proc.time()

vi_agegroup_list <- vector("list", repeat_model)
for (i in 1:repeat_model){
  vi_agegroup_list[[i]] <- list()
}



for (n.repeat_model in 1:repeat_model) {
  # 循环处理每个 explainer_non.fos 对象
  for (n.number_agegroup in 1:dim(Y)[2]) {
    
    # 初始化一个空的数据框，用于累积当前实验的结果
    all_vi_agegroup <- data.frame()
    
    # 这里可以改一下
    #set.seed(i)  # 设置随机种子以确保结果可重复
    vi <- variable_importance(explainer_non.fos_agegroup_list[[n.number_agegroup]][[n.repeat_model]], 
                              loss_function = loss_root_mean_square, B = 100,type="ratio")
    
    # 存储当前年龄组实验的平均变量重要性结果
    vi_agegroup_list[[n.repeat_model]][[n.number_agegroup]] <- vi
    
  }
}

end_time <- proc.time()
elapsed_time <- end_time - start_time
print(elapsed_time)



# 初始化一个列表，用于存储每个年龄组的平均结果
final_mean_vi_agegroup <- vector("list", dim(Y)[2])

# 处理每个年龄组的结果，进行平均聚合
for (n.age_group in 1:dim(Y)[2]) {
  
  age_group_results <- list()
  # age_group_results <- data.frame()  # 初始化为空数据框
  
  # 从每个小列表中提取当前年龄组的数据框并合并
  for (list_idx in 1:repeat_model) {
    age_group_results <- rbind(age_group_results, vi_agegroup_list[[list_idx]][[n.age_group]])
  }
  
  # 按照变量名、置换和标签聚合平均 dropout_loss
  mean_age_group_results <- aggregate(dropout_loss ~ variable + permutation + label, data = age_group_results, FUN = mean)
  class(mean_age_group_results) <- class(vi)
  mean_age_group_results <- mean_age_group_results[, colnames(vi)]
  
  # 按平均 dropout_loss 排序
  mean_age_group_results <- mean_age_group_results[order(mean_age_group_results$dropout_loss, decreasing = TRUE), ]
  
  # 将当前年龄组的平均结果存储到列表中
  final_mean_vi_agegroup[[n.age_group]] <- mean_age_group_results
  
}

#View(final_mean_vi_agegroup)
plot(final_mean_vi_agegroup[[1]])

#save(final_mean_vi_agegroup, file = "保存数据/final_mean_vi_agegroup.RData")


# for(i in 1:dim(Y)[2]){
#   #plot(final_mean_vi_agegroup[[i]])
#   print(final_mean_vi_agegroup[[i]]$variable)
#   cat("\n")
# }

plot(final_mean_vi_agegroup[[1]])
plot(final_mean_vi_agegroup[[2]])
plot(final_mean_vi_agegroup[[3]])
plot(final_mean_vi_agegroup[[4]])
plot(final_mean_vi_agegroup[[5]])
plot(final_mean_vi_agegroup[[6]])
plot(final_mean_vi_agegroup[[7]])
plot(final_mean_vi_agegroup[[8]])
plot(final_mean_vi_agegroup[[9]])
plot(final_mean_vi_agegroup[[10]])
plot(final_mean_vi_agegroup[[11]])
# plot(final_mean_vi_agegroup[[12]])
# plot(final_mean_vi_agegroup[[13]])
# plot(final_mean_vi_agegroup[[14]])
# plot(final_mean_vi_agegroup[[15]])
# plot(final_mean_vi_agegroup[[16]])
# plot(final_mean_vi_agegroup[[17]])
# plot(final_mean_vi_agegroup[[18]])
# plot(final_mean_vi_agegroup[[19]])
# plot(final_mean_vi_agegroup[[20]])
print(final_mean_vi_agegroup[[1]])

# * 整体：分年龄段取平均 -----------------------------------------------------------------

overall_sum <- data.frame()  # 初始化为空数据框

# 从每个小列表中提取当前年龄组的数据框并合并
for (i in 1:dim(Y)[2]) {
  overall_sum <- rbind(overall_sum, final_mean_vi_agegroup[[i]])
}
overall_sum$label <- "model"

# 按照变量名、置换和标签聚合平均 dropout_loss
overall_mean <- aggregate(dropout_loss ~ variable + permutation + label, data = overall_sum, FUN = mean)
class(overall_mean) <- class(vi)
overall_mean <- overall_mean[, colnames(vi)]

# 按平均 dropout_loss 排序
overall_mean <- overall_mean[order(overall_mean$dropout_loss, decreasing = TRUE), ]

# 将当前年龄组的平均结果存储到列表中
plot(overall_mean)
print(overall_mean)

#save(overall_mean, file = "保存数据/overall_mean.RData")


# shap值分析----------------------------------------------------------------

repeat_shap <- 1


# 初始化函数
initialize_list <- function(dimensions) {
  if (length(dimensions) == 1) {
    return(vector("list", length = dimensions[1]))
  }
  return(lapply(1:dimensions[1], function(x) initialize_list(dimensions[-1])))
}

# 初始化 S_list 和 x_list
S_list <- initialize_list(c(number_agegroup, repeat_model, 124))
x_list <- initialize_list(c(number_agegroup, repeat_model, 124))
shapviz_obj <- initialize_list(c(repeat_model, number_agegroup, 124, repeat_shap))


for (n.repeat in 1:repeat_model) {
  
  #set.seed(123)  # 固定随机种子
  
  # 循环处理每个 SHAP 对象
  for (n.agegroup in 1:number_agegroup) {
    
    explainer <- explainer_non.fos_agegroup_list[[n.agegroup]][[n.repeat]]  # 获取当前的 explainer
    
    for (i in 1:124) {
      # 假设 X 是包含观测值的数据框或矩阵，这里用 X[i,] 表示第 i 行观测值
      new_observation <- X[i, ]
      
      # 存储所有 SHAP 值的列表
      all_shap_values_S <- list()
      
      # 循环计算 SHAP 值
      for (j in 1:repeat_shap) {
        shap_values <- predict_parts_shap(explainer, new_observation = new_observation)
        shapviz_obj[[n.repeat]][[n.agegroup]][[i]][[j]] <- 
          shapviz(shap_values)  # 对 shap_values 进行 shapviz
        all_shap_values_S[[j]] <- shapviz_obj[[n.repeat]][[n.agegroup]][[i]][[j]]$S  # 存储 shapviz 后的对象
      }
      
      # 初始化一个和矩阵相同维度的零矩阵
      avg_matrix <- matrix(0, nrow = nrow(all_shap_values_S[[1]]), ncol = ncol(all_shap_values_S[[1]]))
      
      # 逐元素求和
      for (mat in all_shap_values_S) {
        avg_matrix <- avg_matrix + mat
      }
      
      # 求平均
      avg_matrix <- avg_matrix / length(all_shap_values_S)
      
      #  # 使用 sapply 来简化矩阵求和和平均的过程
      # avg_matrix <- do.call(rbind, all_shap_values_S) %>% rowMeans()
      # avg_matrix <- matrix(avg_matrix, nrow = nrow(all_shap_values_S[[1]]), ncol = ncol(all_shap_values_S[[1]]))
      
      # 存储到相应的位置
      S_list[[n.agegroup]][[n.repeat]][[i]] <- avg_matrix
      x_list[[n.agegroup]][[n.repeat]][[i]] <- shapviz_obj[[n.repeat]][[n.agegroup]][[i]][[1]]$X
    }
  }
}


#备份
S_list_0 <- S_list
x_list_0 <- x_list

shapviz_obj_list_baseline  <- vector("list", repeat_model)

for (i in seq_len(repeat_model)) {
  shapviz_obj_list_baseline[[i]] <- vector("list", number_agegroup)
  for (j in seq_len(number_agegroup)) {
    shapviz_obj_list_baseline[[i]][[j]] <- vector("list", 124)
    for (k in seq_len(124)) {
      shapviz_obj_list_baseline[[i]][[j]][[k]] <- vector("list", repeat_shap)
      for (l in seq_len(repeat_shap)) {
        shapviz_obj_list_baseline[[i]][[j]][[k]][[l]] <- shapviz_obj[[i]][[j]][[k]][[l]]$baseline  # 或者初始化为某个值，例如 NA 或一个空列表
      }
    }
  }
}
#View(shapviz_obj_list_baseline)
#不同模型不一样

#View(shapviz_obj_list_baseline[[1]])
#不同shap循环一样
#不同点一样（124）
#不同年龄段不一样
shapviz_obj_list_baseline_0 <- shapviz_obj_list_baseline

#继续处理shapviz_obj_list_baseline，取repeat_model的平均值
shapviz_obj_list_baseline_model <- vector("list", repeat_model)

for (i in 1:repeat_model){
  shapviz_obj_list_baseline_model[[i]] <- vector("list", number_agegroup)
  for (j in 1:number_agegroup){
    shapviz_obj_list_baseline_model[[i]][[j]] <- shapviz_obj_list_baseline[[i]][[j]][[1]][[1]]
  }
  
}
#View(shapviz_obj_list_baseline_model)

# 初始化存放平均值的向量
average_baseline_values <- numeric(number_agegroup)

# 遍历每个年龄组
for (j in seq_len(number_agegroup)) {
  # 初始化用于累积平均值的变量
  sum_values <- 0
  
  # 遍历每个模型
  for (i in seq_len(repeat_model)) {
    # 获取每个子列表中的元素值（唯一的 numeric 值）
    element <- shapviz_obj_list_baseline_model[[i]][[j]][[1]]
    # 累加元素值
    sum_values <- sum_values + element
  }
  
  # 计算平均值并存入结果向量
  average_baseline_values [j] <- sum_values / repeat_model
}

# 查看平均值向量
print(average_baseline_values)



# *分年龄段----

# 现在 S_mean_list 和 X_mean_list 分别包含了每个解释器的平均 S 值和 X 值
# 合并所有的 S 在第二个维度
S_mean_list <- vector("list", length = number_agegroup)

for (agegroup in 1:number_agegroup) {
  # 将每个模型的 S 合并并计算平均值
  combined_S <- Reduce("+", lapply(S_list[[agegroup]], function(S) do.call(rbind, S)))
  S_mean_list[[agegroup]] <- combined_S / repeat_model
}




# 合并所有的 X 在第二个维度
X_mean_list <- vector("list", length = number_agegroup)

for (agegroup in 1:number_agegroup) {
  combined_X_list <- lapply(x_list[[agegroup]], function(X) {
    X_mat <- do.call(rbind, X)
    X_mat <- as.matrix(X_mat)
    X_mat <- apply(X_mat, 2, as.numeric)
    return(X_mat)
  })
  
  combined_X <- Reduce("+", combined_X_list)
  X_mean_list[[agegroup]] <- combined_X / repeat_model
}


# S_mean_list 和 X_mean_list 现在包含了每个年龄组的平均 S 值和 X 值
for (i in 1:dim(Y)[2]){
  # 找到除了 "性别" 和 "城乡" 列之外的所有列的索引
  exclude_cols <- c("性别", "城乡")
  cols_to_convert <- setdiff(names(X_mean_list[[i]]), exclude_cols)
  
  X_mean_list[[i]] <- as.data.frame(X_mean_list[[i]])
  X_mean_list[[i]][, cols_to_convert] <- apply(X_mean_list[[i]][, cols_to_convert], 2, as.numeric)
  
  # 确保 "性别" 和 "城乡" 列保留为字符类型
  X_mean_list[[i]][, exclude_cols] <- lapply(X_mean_list[[i]][, exclude_cols], as.character)
  
}


# 使用 shapviz 函数创建合并后的 shapviz 对象
shapviz_obj_list <- vector("list", length = number_agegroup)

for(i in 1:dim(Y)[2]){
  
  shapviz_obj_list[[i]] <- shapviz(S_mean_list[[i]], X_mean_list[[i]], baseline = average_baseline_values[i])
  
}


# 打印合并后的 shapviz 对象
print(shapviz_obj_list)


# 绘制 SHAP 值的重要性图
sv_importance(shapviz_obj_list[[1]], kind = "beeswarm", show_numbers = TRUE)
#sv_importance(shapviz_obj, show_numbers = TRUE)

sv_dependence(shapviz_obj_list[[1]], v = c("性别"),color_var ="城乡")
sv_importance(shapviz_obj_list[[2]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[3]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[4]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[5]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[6]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[7]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[8]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[9]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[10]], kind = "beeswarm", show_numbers = TRUE)
sv_importance(shapviz_obj_list[[11]], kind = "beeswarm", show_numbers = TRUE)



# *整体模型----

# 使用 Reduce 函数对 20 个矩阵进行元素相加
S_mean_sum <- Reduce(`+`, S_mean_list) / length(S_mean_list)

# 使用 Reduce 函数对 20 个矩阵进行元素相加
X_mean_sum <- Reduce(`+`, X_mean_list) / length(X_mean_list)

#str(X_mean_sum)
# 找到除了 "性别" 和 "城乡" 列之外的所有列的索引
exclude_cols <- c("性别", "城乡")
cols_to_convert <- setdiff(names(X_mean_sum), exclude_cols)

X_mean_sum <- as.data.frame(X_mean_sum)
X_mean_sum[, cols_to_convert] <- apply(X_mean_sum[, cols_to_convert], 2, as.numeric)

# 确保 "性别" 和 "城乡" 列保留为字符类型
X_mean_sum[, exclude_cols] <- lapply(X_mean_sum[, exclude_cols], as.character)

# 使用 shapviz 函数创建合并后的 shapviz 对象
shapviz_obj <- shapviz(S_mean_sum, X_mean_sum,baseline = mean(average_baseline_values))
sv_importance(shapviz_obj, kind = "beeswarm", show_numbers = TRUE)



sv_dependence(shapviz_obj, v = c("性别"),color_var =NULL)
sv_dependence(shapviz_obj, v = c("性别"),color_var ="城乡")
sv_dependence(shapviz_obj, v = c("城乡"),color_var ="性别")

sv_dependence(shapviz_obj, v = vec[3],color_var =c("性别","城乡"))
sv_dependence(shapviz_obj, v = vec[4],color_var =c("性别","城乡"))
sv_dependence(shapviz_obj, v = vec[5],color_var =c("性别","城乡"))
sv_dependence(shapviz_obj, v = vec[6],color_var =c("性别","城乡"))
sv_dependence(shapviz_obj, v = vec[7],color_var =c("性别","城乡"))
sv_dependence(shapviz_obj, v = vec[8],color_var =c("性别","城乡"))

sv_dependence(shapviz_obj, v = c("性别","城乡"),color_var =vec[3])
sv_dependence(shapviz_obj, v = c("性别","城乡"),color_var =vec[4])
sv_dependence(shapviz_obj, v = c("性别","城乡"),color_var =vec[5])
sv_dependence(shapviz_obj, v = c("性别","城乡"),color_var =vec[6])
sv_dependence(shapviz_obj, v = c("性别","城乡"),color_var =vec[7])
sv_dependence(shapviz_obj, v = c("性别","城乡"),color_var =vec[8])




# 曲线图分析----------------------------------------------------------------

# *PDP----
# 初始化存储部分依赖图结果的列表
pdp_agegroup_list <- list()

# 外层循环：对每个年龄组进行操作
for (i in 1:number_agegroup) {
  
  # 创建存储该年龄组每次建模结果的列表
  pdp_agegroup_list[[i]] <- vector("list", repeat_model)
  
  # 对每次建模计算部分依赖图
  for (n.repeat_model in 1:repeat_model) {
    explainer_non.fos_agegroup <- explainer_non.fos_agegroup_list[[i]][[n.repeat_model]]
    
    pdp_agegroup <- variable_effect(explainer_non.fos_agegroup,
                                    variables = colnames(X)[8],
                                    type = "partial_dependency")
    
    # 将部分依赖图结果存储在对应的位置
    pdp_agegroup_list[[i]][[n.repeat_model]] <- pdp_agegroup
  }
}

# pdp_agegroup_list[[1]][[1]] ==pdp_agegroup_list[[1]][[2]]


pdp_avg <- vector("list", length = number_agegroup)

# 将每个数据框的 _yhat_ 列填充到 yhat_matrix 中
for (i in 1:number_agegroup){
  
  yhat_matrix <- matrix(0, nrow = nrow(pdp_agegroup_list[[i]][[1]]), ncol = repeat_model)
  
  for (j in 1:repeat_model) {
    
    yhat_matrix[, j] <- pdp_agegroup_list[[i]][[j]]$`_yhat_`
    
  }
  # 计算每行 _yhat_ 的均值
  yhat_avg <- rowMeans(yhat_matrix)
  
  # 创建新的数据框，并将 _yhat_ 列替换为均值
  pdp_avg[[i]] <- pdp_agegroup_list[[i]][[1]]
  pdp_avg[[i]]$`_yhat_` <- yhat_avg
  
}

plot(pdp_avg[[1]])


# *ALE----

ale_agegroup_list <- list()

# 外层循环：对每个年龄组进行操作
for (i in 1:number_agegroup) {
  
  # 创建存储该年龄组每次建模结果的列表
  ale_agegroup_list[[i]] <- vector("list", repeat_model)
  
  for (n.repeat_model in 1:repeat_model) {
    explainer_non.fos_agegroup <- explainer_non.fos_agegroup_list[[i]][[n.repeat_model]]
    
    ale_agegroup <- variable_effect(explainer_non.fos_agegroup,
                                    variables = colnames(X)[8],
                                    type = "accumulated_dependency")
    
    ale_agegroup_list[[i]][[n.repeat_model]] <- ale_agegroup
  }
}

#ale_agegroup_list[[1]][[1]] ==ale_agegroup_list[[1]][[2]]


ale_avg <- vector("list", length = number_agegroup)

# 将每个数据框的 _yhat_ 列填充到 yhat_matrix 中
for (i in 1:number_agegroup){
  
  yhat_matrix <- matrix(0, nrow = nrow(ale_agegroup_list[[i]][[1]]), ncol = repeat_model)
  
  for (j in 1:repeat_model) {
    
    yhat_matrix[, j] <- ale_agegroup_list[[i]][[j]]$`_yhat_`
    
  }
  # 计算每行 _yhat_ 的均值
  yhat_avg <- rowMeans(yhat_matrix)
  
  # 创建新的数据框，并将 _yhat_ 列替换为均值
  ale_avg[[i]] <- ale_agegroup_list[[i]][[1]]
  ale_avg[[i]]$`_yhat_` <- yhat_avg
  
}

plot(ale_avg[[1]])
plot(ale_avg)


