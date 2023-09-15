# TO RUN: download the data folder and change the location in the setwd() command to the parent directory of the data folder
# Windows users will need to change "data/" to "data\\" throughout the script

packages <- c("tidyverse")
install.packages(setdiff(packages, rownames(installed.packages())))

library(tidyverse)

# change the working directory as needed. The working directory should be the parent folder of the 'data' folder. 
# please note that paths require double backslashes when running scripts on Windows.
setwd("/Users/xyz")

b2AR_Gs_data <- read_csv("data/b2AR_Gs_data.csv")
residue_classification <- read_csv("data/residue_classification.csv")
R_ligand_4LDO <- read_tsv("data/4LDO_A.txt")
R_NAM <- read_tsv("data/6OBA_A.txt")
R_PAM <- read_tsv("data/6N48_A.txt")


# adrenaline/epinephrine binding

R_ligand_contacts_4LDO_intermediate <- R_ligand_4LDO %>%
  filter(Res1 == "ALE") %>%
  rename(atomic_contacts = `Number of atomic contacts`) %>%
  mutate(ResNum2 = ResNum2-1000) %>%
  separate(Atoms, into = c("Atom1", "Atom2"),sep ="-") 


# ligand contacts for Figure 3B - details
R_ligand_contacts_4LDO <- R_ligand_contacts_4LDO_intermediate%>%
  group_by(Res1, ResNum1, Chain1, Res2, ResNum2, Atom1) %>%
  summarise(contacts = n()) %>%
  left_join(R_ligand_contacts_4LDO_intermediate %>% select(ResNum2, Atom1, Atom2), by = c("ResNum2", "Atom1")) %>%
  left_join(b2AR_Gs_data %>% 
              select(amino_acid, GPCRdb, efficacy, potency, expression), by = c("ResNum2" = "amino_acid")) %>%
  mutate(R_position = str_c(Res2, ResNum2), PDB = "4LDO") %>%
  select(PDB, everything())

# ligand contacts for Figure 3B - summary
Figure_3B <- R_ligand_contacts_4LDO %>%
  # to get receptor residue - atom contact numbers only 
  group_by(ResNum2, Atom1, contacts) %>%
  summarise()


# NAM binding

R_NAM_contacts_intermediate <- R_NAM %>%
  filter(Res1 == "M3J") %>%
  rename(atomic_contacts = `Number of atomic contacts`) %>%
  separate(Atoms, into = c("Atom1", "Atom2"),sep ="-") 

R_NAM_contacts <- R_NAM_contacts_intermediate%>%
  group_by(Res1, ResNum1, Chain1, Res2, ResNum2, Atom1) %>%
  summarise(contacts = n()) %>%
  left_join(R_NAM_contacts_intermediate %>% select(ResNum2, Atom1, Atom2), by = c("ResNum2", "Atom1")) %>%
  left_join(b2AR_Gs_data %>% 
              select(amino_acid, GPCRdb, potency, efficacy, expression), by = c("ResNum2" = "amino_acid")) %>%
  left_join(residue_classification %>% select(amino_acid, classification), by = c("ResNum2" = "amino_acid")) %>%
  mutate(R_position = str_c(Res2, ResNum2)) %>%
  mutate(Atom1 = factor(Atom1, levels = c("C15", "C16", "C14", "C18", "C13", "C19", "N12","C02", "N01", "C04", "N03", "N05", "C06", "C11", "C07", "C10", "C08", "C09"))) %>%
  arrange(ResNum2)

Figure_6E <- R_NAM_contacts %>%
  select(-(Atom2)) %>%
  unique()

# PAM binding

R_PAM_contacts_intermediate <- R_PAM %>%
  filter(Res1 == "KBY") %>%
  rename(atomic_contacts = `Number of atomic contacts`) %>%
  mutate(ResNum2 = ResNum2-1000) %>%
  separate(Atoms, into = c("Atom1", "Atom2"),sep ="-") 

R_PAM_contacts <- R_PAM_contacts_intermediate%>%
  group_by(Res1, ResNum1, Chain1, Res2, ResNum2, Atom1) %>%
  summarise(contacts = n()) %>%
  ungroup()%>%
  left_join(R_PAM_contacts_intermediate %>% select(ResNum2, Atom1, Atom2), by = c("ResNum2", "Atom1")) %>%
  left_join(b2AR_Gs_data %>% 
              select(amino_acid, GPCRdb, potency, efficacy, expression), by = c("ResNum2" = "amino_acid")) %>%
  mutate(R_position = str_c(Res2, ResNum2)) %>%
  left_join(residue_classification %>% select(amino_acid, classification), by = c("ResNum2" = "amino_acid")) %>%
  mutate(Atom1 = factor(Atom1, levels = c("C16", "O3", "C2", "C14", "C11", "C10", "C4", "C3","C17", "O4", "C13", "C19", "C12", "C1", "C6", "C8", "C9", "C30", "C31", "C26", "C27")))%>%
  arrange(ResNum2)

Figure_6F <- R_PAM_contacts %>%
  # to get receptor residue - atom contact numbers only 
  select(-(Atom2)) %>%
  unique()