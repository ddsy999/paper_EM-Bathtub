colnames(theta_df) <- column_names

# theta_df
p1=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=sumQfunc %>% as.numeric()))+geom_point()+geom_line()+ggtitle("sumQfunc")
p2=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=beta1 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Beta1")
p3=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=beta3 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Beta3")

p4=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=diffB_beta1 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("diffB_beta1")
p5=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=diffB_beta3 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("diffB_beta3")

p6=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=Q1 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Q1")
p7=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=Q2 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Q2")
p8=theta_df %>% ggplot(aes(x=iter%>% as.numeric(),y=Q3 %>% as.numeric()))+geom_point()+geom_line()+ggtitle("Q3")
# 두 개의 그림을 결합하여 하나로 표시
combined_plot <- (p1+p2+p3)/(p4+p5)/(p6+p7+p8)

# 결합된 그림을 출력
print(combined_plot)



initial_latent_varialbe = result_latentZ_mat[[1]]
middle_latent_varialbe = result_latentZ_mat[[floor(length(result_latentZ_mat)/2)]]
final_latent_varialbe = result_latentZ_mat[[length(result_latentZ_mat)]]

colnames(initial_latent_varialbe)=c("Z1","Z2","Z3")
colnames(middle_latent_varialbe)=c("Z1","Z2","Z3")
colnames(final_latent_varialbe)=c("Z1","Z2","Z3")

ggplot(middle_latent_varialbe,aes(x=time_vec,y=Z1),color="black")+geom_line()+
geom_line(data=middle_latent_varialbe,aes(x=time_vec,y=Z3),color="black")+
geom_line(data=final_latent_varialbe,aes(x=time_vec,y=Z1),color="red")+
geom_line(data=final_latent_varialbe,aes(x=time_vec,y=Z3),color="red")










