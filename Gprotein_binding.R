# TO RUN: download the data folder and change the location in the setwd() command to the parent directory of the data folder
# Windows users will need to change "data/" to "data\\" throughout the script
# Common G protein Numbering (CGN) was obtained from the CGN server, https://www.mrc-lmb.cam.ac.uk/CGN/

packages <- c("tidyverse")
install.packages(setdiff(packages, rownames(installed.packages())))

library(tidyverse)

# change the working directory as needed. The working directory should be the parent folder of the 'data' folder. 
# please note that paths require double backslashes when running scripts on Windows.
setwd("/Users/xyz")

b2AR_Gs_data <- read_csv("data/b2AR_Gs_data.csv")
R_alpha_contacts <- read_tsv("data/3SN6_A-R.txt")
R_beta_contacts <- read_tsv("data/3SN6_B-R.txt")

CGN <- read_tsv("data/3SN6_CGN.txt") %>% 
  separate(`3SN6`, into = c("Res1", "ResNum1"),sep = 1) %>% 
  select(-`Sort number`)%>%
  mutate(ResNum1 = as.numeric(ResNum1))

R_Gs_contacts <- R_alpha_contacts%>%
  bind_rows(R_beta_contacts) 

# compare with Figure 3E
R_Gs_contacts_summary <- R_alpha_contacts%>%
  bind_rows(R_beta_contacts) %>%
  rename(atomic_contacts = `Number of atomic contacts`) %>%
  group_by(Res1, ResNum1, Chain1, Res2, ResNum2, atomic_contacts) %>%
  summarise() %>%
  left_join(CGN, by = c("Res1", "ResNum1")) %>%
  left_join(b2AR_Gs_data %>% 
              select(amino_acid, GPCRdb, efficacy, potency), by = c("ResNum2" = "amino_acid")) %>%
  replace_na(list(CGN = "G beta")) %>%
  mutate(G_position = str_c(Res1, ResNum1), R_position = str_c(Res2, ResNum2)) %>%
  arrange(ResNum2, ResNum1)