

explained_variance_numeric2(maf=my_data_harm$eaf.bmi, 
                            beta=my_data_harm$beta.bmi, 
                            se_beta=my_data_harm$se.bmi, 
                            samplesize=n[length(n)]) %>% sum


explained_variance_numeric2(maf=my_data_harm$eaf.bmi, 
                            beta=my_data_harm$beta.bmi, 
                            se_beta=my_data_harm$se.bmi, 
                            samplesize=1000) %>% sum


explained_variance_numeric2(maf=my_data_harm$eaf.bmi[my_data_harm$p_value_exp_sample_100000 < 5*10^(-6)], 
                            beta=my_data_harm$beta.bmi[my_data_harm$p_value_exp_sample_100000 < 5*10^(-6)], 
                            se_beta=my_data_harm$se.bmi[my_data_harm$p_value_exp_sample_100000 < 5*10^(-6)], 
                            samplesize=100000) %>% sum


explained_variance_numeric2(maf=my_data_harm$eaf.bmi[my_data_harm$p_value_exp_sample_150000 < 5*10^(-6)], 
                            beta=my_data_harm$beta.bmi[my_data_harm$p_value_exp_sample_150000 < 5*10^(-6)], 
                            se_beta=my_data_harm$se.bmi[my_data_harm$p_value_exp_sample_150000 < 5*10^(-6)], 
                            samplesize=150000) %>% sum


