
# PLAUSIBILITY CHECK AFTER DATA HARMONIZATION

my_data %>% 
  left_join(my_data_harm, by=c("SNP"="SNP") ) -> test

test %>% select(beta_exp, beta.bmi, se_exp, se.bmi) -> test1
sum(test1$beta_exp == test1$beta.bmi, na.rm = TRUE)
sum(test1$se_exp == test1$se.bmi, na.rm = TRUE)
## IDENTICAL

test %>% 
  left_join(GBC_dat,by=c("SNP"="SNP")) -> test2

test2 %>% select(beta_out.y, beta.gbc, se_out.y, se.gbc) -> test3
sum(test3$beta_out.y == test3$beta.gbc, na.rm = TRUE)
sum(test3$se_out.y == test3$se.gbc, na.rm = TRUE)
