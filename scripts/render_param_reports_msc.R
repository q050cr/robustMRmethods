



## render report and mail results
library(dplyr)
library(gmailr)
rm(list=ls())  # not in rmd doc, otherwise params deleted
#gm_auth_configure(path = "/Users/christophreich/Library/CloudStorage/OneDrive-bwedu/credentials.json")
#gm_auth()

render_report = function(...) {
  ## disease: acs, cad, hfref, dcm | default: dcm
  ## analysis: selected, full | default: selected
  
  rmarkdown::render(
    "scripts/MSc.Rmd",
    params = list(
      doc_title =  paste0("Assessment of simulation procedures and robust methods\nfor Mendelian randomization analysis")
    ),
    output_file = paste0(
      format(Sys.time(), "%Y%m%d"), '_robustMRsimulation.html'
    ),
    output_dir = "reports",
    envir = globalenv()
  )
  # save filename for sending later
  filename.html <<-  paste0("reports/", format(Sys.time(), "%Y%m%d"), '_robustMRsimulation.html' )
  title.mail <<- paste0("Robust MR report | Date:  ",  format(Sys.time(), "%Y-%m-%d"))
  
 # email <-
 #   gmailr::gm_mime() %>%
 #   gmailr::gm_to("reich.c@hotmail.com") %>%
 #   gmailr::gm_from("q050cr@gmail.com") %>%
 #   gmailr::gm_subject(paste0("Job finished: ", title.mail)) %>%
 #   gmailr::gm_text_body("Hi Christoph,\nHere is a new report for you.\nCheers")
 # # attachment
 # email <- gmailr::gm_attach_file(email, filename.html)
 # gmailr::gm_send_message(email, user_id = "me")
}

# NOW START!!!
render_report()
