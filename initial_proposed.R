source("DAEM_BarrierMethod_function.R")

# Data 설정
file_name = 'FRT_censord.txt'
fdata = read.table(file_name, header = TRUE)
dataName = tools::file_path_sans_ext(file_name)

N = nrow(fdata)
k = 3
event_vec = as.numeric(fdata[,2])
time_vec = as.numeric(fdata[,1])

tot = 1e-9
maxEMIter = 1e+5
maxIterAnnealing = 100
learningRateBp = 2

annealingSchedule = seq(0.2, 0.999999999, length.out = maxIterAnnealing)
annealingSchedule = rep(1,maxIterAnnealing)
bpBaseSchedule = exp(seq(log(1e+2), log(1e+7), length.out = maxIterAnnealing))

# 초기값 리스트 정의
initial_beta_list = list(
  c(0.1, 1, 1.5),
  c(0.5, 1 ,2),
  c(0.9, 1, 10)
)

initial_pi_list = list(
  c(1, 1, 1),
  c(1, 2, 1),
  c(1, 8, 1)
)

# 반복 수행
for (b_id in seq_along(initial_beta_list)) {
  for (p_id in seq_along(initial_pi_list)) {
    
    cat("\n===== Run with beta", b_id, "and pi", p_id, "=====\n")
    
    theta_df_full = NULL
    result_latentZ_mat = list()
    
    initial_beta = initial_beta_list[[b_id]]
    initial_pi_set = initial_pi_list[[p_id]]
    initial_pi = initial_pi_set / sum(initial_pi_set)
    
    initial_lambda = initial_lambda_func(time_vec, event_vec,
                                         initial_beta, ratio1 = 0.3, ratio3 = 0.7)
    
    beta_vec = initial_beta
    pi_vec = initial_pi
    lambda_vec = initial_lambda
    latentZ_mat = Estep_result(beta_vec, lambda_vec, pi_vec, alpha = 1)
    
    for (ITerAnneal in 1:maxIterAnnealing) {
      annealingPara = annealingSchedule[ITerAnneal]
      bpBase = bpBaseSchedule[ITerAnneal]
      
      for (iter in 1:maxEMIter) {
        new_pi = colSums(latentZ_mat) / N
        
        candi_before_vec = beta_vec
        bp1 = bpBase
        bp3 = bpBase
        
        new_beta1 = barrier_beta1(candi_before_vec[1], latentZ_mat, bp = bpBase)
        new_beta3 = barrier_beta3(candi_before_vec[3], latentZ_mat, bp = bpBase)
        new_beta = c(new_beta1, 1, new_beta3)
        
        new_lambda = sapply(1:k, function(i)
          sum(latentZ_mat[, i] * event_vec) / sum(latentZ_mat[, i] * (time_vec ^ new_beta[i]))
        )
        
        parameter_diff = sqrt(sum((beta_vec - new_beta)^2 + (pi_vec - new_pi)^2))
        
        beta_vec = new_beta
        pi_vec = new_pi
        lambda_vec = new_lambda
        alpha_temper = annealingPara
        
        if (parameter_diff < tot || iter == maxEMIter) {
          printResult()
          print("!!! parameter diff Break !!!")
          
          theta_df_full = rbind(theta_df_full,
                                data.frame(
                                  data_Name = dataName,
                                  Iter = ITerAnneal,
                                  beta1 = beta_vec[1],
                                  beta3 = beta_vec[3],
                                  lambda1 = lambda_vec[1],
                                  lambda2 = lambda_vec[2],
                                  lambda3 = lambda_vec[3],
                                  diffBeta1 = diffB_onlyB(beta_vec[1], latentZ_mat, j = 1),
                                  diffBeta3 = diffB_onlyB(beta_vec[3], latentZ_mat, j = 3),
                                  diffLambda1 = diffL(beta_vec[1], lambda_vec[1], latentZ_mat, j = 1),
                                  diffLambda3 = diffL(beta_vec[3], lambda_vec[3], latentZ_mat, j = 3),
                                  pi1 = pi_vec[1],
                                  pi2 = pi_vec[2],
                                  pi3 = pi_vec[3],
                                  init_beta1 = initial_beta[1],
                                  init_beta2 = initial_beta[2],
                                  init_beta3 = initial_beta[3],
                                  init_pi1 = initial_pi[1],
                                  init_pi2 = initial_pi[2],
                                  init_pi3 = initial_pi[3],
                                  EM_endIter = iter,
                                  tempering = alpha_temper,
                                  Qlike = sumQfunc(beta_vec, lambda_vec, latentZ_mat),
                                  Q1 = Qfunc_onlyB(beta = beta_vec[1], latentZ_mat, j = 1),
                                  Q2 = Qfunc_onlyB(beta = beta_vec[2], latentZ_mat, j = 2),
                                  Q3 = Qfunc_onlyB(beta = beta_vec[3], latentZ_mat, j = 3),
                                  bpBase = bpBase,
                                  bp1 = bp1,
                                  bp3 = bp3
                                )
          )
          
          result_latentZ_mat[[ITerAnneal]] = latentZ_mat
          break
        }
        
        latentZ_mat = Estep_result(beta_vec, lambda_vec, pi_vec, alpha = alpha_temper)
      }
    }
    
    # 결과 저장
    file_tag = paste0("b", b_id, "_p", p_id)
    # result_Name = paste0("ProposedResult_", dataName, "_", file_tag, ".txt")
    result_Name = paste0("ProposedResult_noDAEM_", dataName, "_", file_tag, ".txt")
    write.table(theta_df_full, file = result_Name, row.names = FALSE)
    
    cat("Saved to", result_Name, "\n")
  }
}








# 모든 txt 파일 불러오기
file_pattern <- paste0("^ProposedResult_", dataName, "_.*\\.txt$")
file_list <- list.files(pattern = file_pattern)

# 결과 저장할 데이터프레임
result_summary <- data.frame()

# 각 파일에 대해 반복
for (file_name in file_list) {
  
  # 데이터 불러오기
  df <- read.table(file_name, header = TRUE)
  
  # abs(diffBeta1) + abs(diffBeta3) 계산
  df <- df %>%
    mutate(abs_diff_sum = abs(diffBeta1) + abs(diffBeta3))
  
  # 최소값 행 추출
  min_row <- df %>%
    filter(abs_diff_sum == min(abs_diff_sum)) %>%
    mutate(source_file = file_name)
  
  # 결과 누적
  result_summary <- bind_rows(result_summary, min_row)
}

# 정렬된 결과 확인
result_summary_sorted <- result_summary %>%
  select(source_file, Iter, beta1, beta3, diffBeta1, diffBeta3, abs_diff_sum) %>%
  arrange(abs_diff_sum)

# 출력
print(result_summary_sorted)










# 모든 txt 파일 불러오기
file_list <- list.files(pattern = "^ProposedResult_noDAEM_.*\\.txt$")

# 결과 저장할 데이터프레임
result_summary <- data.frame()

# 각 파일에 대해 반복
for (file_name in file_list) {
  
  # 데이터 불러오기
  df <- read.table(file_name, header = TRUE)
  
  # abs(diffBeta1) + abs(diffBeta3) 계산
  df <- df %>%
    mutate(abs_diff_sum = abs(diffBeta1) + abs(diffBeta3))
  
  # 최소값 행 추출
  min_row <- df %>%
    filter(abs_diff_sum == min(abs_diff_sum)) %>%
    mutate(source_file = file_name)
  
  # 결과 누적
  result_summary <- bind_rows(result_summary, min_row)
}

# 정렬된 결과 확인
result_summary_sorted <- result_summary %>%
  select(source_file, Iter, beta1, beta3, diffBeta1, diffBeta3, abs_diff_sum) %>%
  arrange(abs_diff_sum)

# 출력
print(result_summary_sorted)


read.table("ProposedResult_noDAEM_FRT_censord_b2_p2.txt", header = TRUE) %>% pull(diffBeta1) %>% plot
read.table("ProposedResult_noDAEM_FRT_censord_b2_p2.txt", header = TRUE) %>% pull(diffBeta3) %>% plot
read.table("ProposedResult_noDAEM_FRT_censord_b2_p2.txt", header = TRUE) %>% pull(Qlike) %>% plot


read.table("ProposedResult_noDAEM_FRT_censord_b3_p1.txt", header = TRUE) %>% pull(diffBeta1) %>% plot
read.table("ProposedResult_noDAEM_FRT_censord_b3_p2.txt", header = TRUE) %>% pull(diffBeta3) %>% plot
read.table("ProposedResult_noDAEM_FRT_censord_b2_p3.txt", header = TRUE) %>% pull(Qlike) %>% plot


