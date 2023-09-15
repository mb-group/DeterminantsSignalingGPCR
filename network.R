# The required .txt files were downloaded from https://pca.mbgroup.bio/index.html using default settings
# The analysis is based on 2RH1 - an inactive-state structure of b2AR and 3SN6 - an active-state structure of b2AR in complex with Gs
# TO RUN: download the data folder and change the location in the setwd() command to the parent directory of the data folder
# Windows users will need to change "data/" to "data\\" in the read_tsv() and read_xlsx() function calls  

# install required packages
packages <- c("tidyverse")
install.packages(setdiff(packages, rownames(installed.packages())))

library(tidyverse)

# change the working directory as needed. The working directory should be the parent folder of the 'data' folder. 
# please note that paths require double backslashes when running scripts on Windows.

setwd("/Users/xyz")


# import contacts from Protein Contacts Atlas
contacts_2RH1 <- read_tsv("data/2RH1_A.txt")
contacts_3SN6 <- read_tsv("data/3SN6_R.txt")

# import table with signalling data
b2AR_Gs_data <- read_csv("data/b2AR_Gs_data.csv")


# function for cleaning and annotating contacts .txt files downloaded from Protein Contacts Atlas
# note 1: this function needs to be adapted to work with PDBs where the amino acid number and the residue number do not match
# note 2: this function needs to be adapted for structures with more than 399 resolved GPCR residues
clean_contacts <- as_mapper(~.x %>%
      filter(ResNum1 < 400 & ResNum2 < 400) %>%
      left_join(b2AR_Gs_data %>% select(amino_acid, SSE, GPCRdb) %>% rename(ResNum1 = amino_acid, SSE_1 = SSE, GPCRdb_1 = GPCRdb), by = "ResNum1") %>%
      left_join(b2AR_Gs_data %>% select(amino_acid, SSE, GPCRdb) %>% rename(ResNum2 = amino_acid, SSE_2 = SSE, GPCRdb_2 = GPCRdb), by = "ResNum2") %>%
      filter(`Chain Types` != "M-M") %>%
      filter(SSE_1 != SSE_2) %>%
      separate(Atoms, into = c("Atom1", "Atom2"), sep = "-", remove = FALSE))

# clean contacts of 2RH1 and 3SN6
contacts_2RH1_clean <- clean_contacts(contacts_2RH1)
contacts_3SN6_clean <- clean_contacts(contacts_3SN6)

# residues resolved in one or both structures (2RH1 & 3SN6)
resolved_residues <- tibble(amino_acid = 2:413) %>%
  mutate(resolved_in_2RH1 = case_when(
    amino_acid %in% c(contacts_2RH1$ResNum1[contacts_2RH1$ResNum1 %in% c(1:342)], 
                      contacts_2RH1$ResNum2[contacts_2RH1$ResNum2 %in% c(1:342)]) ~1,
    TRUE ~0
  ), resolved_in_3SN6 = case_when(
    amino_acid %in% c(contacts_3SN6$ResNum1[contacts_3SN6$ResNum1 %in% c(1:342)], 
                      contacts_3SN6$ResNum2[contacts_3SN6$ResNum2 %in% c(1:342)]) ~1,
    TRUE ~0
  )) %>%
  left_join(b2AR_Gs_data %>% select (amino_acid, GPCRdb), by = "amino_acid")

# overview (count) of residues resolved in the two structures
resolved_residues %>%
  group_by(resolved_in_2RH1, resolved_in_3SN6) %>%
  count()

# residues resolved in both structures (2RH1 & 3SN6)
resolved_in_2RH1_3SN6 <- resolved_residues %>%
  filter(resolved_in_2RH1 == 1 & resolved_in_3SN6 == 1)

# combining contacts of 2RH1 and 3SN6
contacts <- contacts_3SN6_clean %>% 
  full_join(contacts_2RH1_clean) %>% 
  select(-(Chain1:SS2), -Distance, -PDB) %>%
  distinct() %>%
  left_join(contacts_2RH1_clean %>% select(ResNum1, ResNum2, Atoms, Distance), by = c("ResNum1", "ResNum2", "Atoms")) %>%
  rename(distance_2RH1 = Distance) %>%
  left_join(contacts_3SN6_clean %>% select(ResNum1, ResNum2, Atoms, Distance), by = c("ResNum1", "ResNum2", "Atoms")) %>%
  rename(distance_3SN6 = Distance) %>%
  mutate(in_2RH1 = case_when(
    distance_2RH1 > 0 ~ 1,
    TRUE ~ 0
  ), in_3SN6 = case_when(
    distance_3SN6 > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  # added to limit the analysis to those residues that are resolved in both 2RH1 and 3SN6
  filter(ResNum1 %in% c(resolved_in_2RH1_3SN6$amino_acid) & ResNum1 %in% c(resolved_in_2RH1_3SN6$amino_acid))

# results for Figure 1C: state-specific and shared contacts
contact_count_Figure1C <- contacts %>%
  group_by(ResNum1, ResNum2)%>%
  summarise(contacts = n(), sum_2RH1 = sum(in_2RH1), sum_3SN6 = sum(in_3SN6)) %>%
  mutate(is_in_2RH1 = case_when(
    sum_2RH1 > 0 ~ 1,
    TRUE ~0
  ), is_in_3SN6 = case_when(
    sum_3SN6 > 0 ~ 1,
    TRUE ~0
  )) %>%
  ungroup() %>%
  mutate(only_in = case_when(
    (is_in_2RH1 == 1 & is_in_3SN6 == 0) ~ "2RH1",
    (is_in_2RH1 == 0 & is_in_3SN6 == 1) ~ "3SN6",
    TRUE ~ "both"
  )) %>%
  filter(ResNum1 %in% c(resolved_in_2RH1_3SN6$amino_acid) & ResNum2 %in% c(resolved_in_2RH1_3SN6$amino_acid))

# mutations: 37 low potency, 21 low efficacy, 21 low potency & efficacy, 3 without measurable signalling
# compare main text, subheading 'One-fifth of the positions in β2AR are important for adrenaline efficacy, potency, or both' 

# get all residues with reduced potency. 
# Total 61 residues: 37 low potency, 21 low potency & efficacy, 3 without measurable signalling (see curation == 2)

Gs_potency_reduced  <- b2AR_Gs_data %>%
  filter(potency == "reduced") %>%
  filter(expression > 25)

# compare with Figures 2A and 2B.
# NB. curation == 2 cannot be plotted since signalling was not measurable. See 2B "no signaling" in this case
Gs_potency_reduced %>%
  group_by(potency, efficacy, curation) %>%
  count()
  

# get all residues with reduced efficacy
# Total 45 residues: 21 low efficacy, 21 low potency & efficacy, 3 without measurable signalling
Gs_efficacy_reduced  <- b2AR_Gs_data %>%
  filter(efficacy == "reduced") %>%
  filter(expression > 25)

# compare with Figures 2A and 2B.
# NB. curation == 2 cannot be plotted since signalling was not measurable. See 2B "no signaling" in this case
Gs_efficacy_reduced %>%  
  group_by(potency, efficacy, curation) %>%
  count()



# get contacts between residues that negatively affect potency upon mutation
contacts_Gs_potency_reduced <- contacts %>%
  filter(ResNum1 %in% Gs_potency_reduced$amino_acid & ResNum2 %in% Gs_potency_reduced$amino_acid)%>%
  group_by(ResNum1, ResNum2, GPCRdb_1, GPCRdb_2) %>%
  summarise(contacts = n(), sum_2RH1 = sum(in_2RH1), sum_3SN6 = sum(in_3SN6)) %>%
  select(GPCRdb_1, GPCRdb_2, everything()) %>%
  mutate(is_in_2RH1 = case_when(
    sum_2RH1 > 0 ~ 1,
    TRUE ~0
  ), is_in_3SN6 = case_when(
    sum_3SN6 > 0 ~ 1,
    TRUE ~0
  ), only_in = case_when(
    (is_in_2RH1 == 1 & is_in_3SN6 == 0) ~ "2RH1",
    (is_in_2RH1 == 0 & is_in_3SN6 == 1) ~ "3SN6",
    TRUE ~ "both"
  ))

# get contacts between residues that negatively affect efficacy upon mutation
contacts_Gs_efficacy_reduced <- contacts %>%
  filter(ResNum1 %in% Gs_efficacy_reduced$amino_acid & ResNum2 %in% Gs_efficacy_reduced$amino_acid)%>%
  group_by(ResNum1, ResNum2, GPCRdb_1, GPCRdb_2) %>%
  summarise(contacts = n(), sum_2RH1 = sum(in_2RH1), sum_3SN6 = sum(in_3SN6)) %>%
  select(GPCRdb_1, GPCRdb_2, everything()) %>%
  mutate(is_in_2RH1 = case_when(
    sum_2RH1 > 0 ~ 1,
    TRUE ~0
  ), is_in_3SN6 = case_when(
    sum_3SN6 > 0 ~ 1,
    TRUE ~0
  ), only_in = case_when(
    (is_in_2RH1 == 1 & is_in_3SN6 == 0) ~ "2RH1",
    (is_in_2RH1 == 0 & is_in_3SN6 == 1) ~ "3SN6",
    TRUE ~ "both"
  ))

# combine contacts between residues negatively affecting potency with those negatively affecting efficacy upon mutation
network <- contacts_Gs_potency_reduced %>%
  bind_rows(contacts_Gs_efficacy_reduced)%>%
  unique() %>%
  ungroup()

# residues and their contacts from Figure 5C/D
# number of active-state specific contacts in the network
network_contacts <- network %>% 
  filter(only_in == "3SN6")

# number of contacts in the network (active-state specific contacts)
# compare with figure 5C
network_contacts %>%
  count()

# number of residues in the network
network_residues <- network %>% 
  filter(only_in == "3SN6") %>%
  select(amino_acid = ResNum1, GPCRdb = GPCRdb_1) %>%
  bind_rows(network %>% 
              filter(only_in == "3SN6") %>%
              select(amino_acid = ResNum2, GPCRdb = GPCRdb_2)) %>%
  unique()

# number of residues in the network
# compare with figure 5C
network_residues %>%
  count()


contacts_summarised <- contact_count_Figure1C %>%
  left_join(b2AR_Gs_data%>% select(amino_acid, GPCRdb_1 = GPCRdb), by = c("ResNum1" = "amino_acid")) %>%
  left_join(b2AR_Gs_data %>% select(amino_acid, GPCRdb_2 = GPCRdb), by = c("ResNum2" = "amino_acid")) %>%
  select(GPCRdb_1, GPCRdb_2, everything())


# all active state-specific contacts
  active_contacts <- contacts_summarised %>%
  filter(only_in == "3SN6") 
  

forms_active_contacts <- tibble(amino_acid = 2:413) %>% 
  left_join(b2AR_Gs_data %>% select(amino_acid, GPCRdb), by = "amino_acid") %>%
  mutate(forms_active_contact = case_when(
    amino_acid %in% c(active_contacts$ResNum1, active_contacts$ResNum2) ~ "yes",
    TRUE ~ "no"
  )) %>%
  mutate(forms_active_contact = factor(forms_active_contact))%>%
  mutate(in_network = case_when(
    amino_acid %in% c(network_contacts$ResNum1, network_contacts$ResNum2) ~ "yes",
    TRUE ~ "no"
  )) %>%
  mutate(active_contact_in_network = case_when(
    (in_network == "yes") ~ "yes, in network",
    (forms_active_contact == "yes") ~"yes, not in network",
    TRUE ~ "no"
  )) %>%
  mutate(active_contact_in_network = factor(active_contact_in_network, levels = c("no", "yes, not in network", "yes, in network")))


residue_classification <- forms_active_contacts %>%
  left_join(b2AR_Gs_data %>% select(pharma_important, amino_acid, expression), by = "amino_acid")%>%
  mutate(classification = case_when(
    expression < 25  ~ NA_character_,
    in_network == "yes" ~ "connected driver",
    (forms_active_contact == "yes" & pharma_important == 1) ~ "disconnected driver",
    (forms_active_contact == "no" & pharma_important == 1) ~ "modulator",
    (forms_active_contact == "yes" & pharma_important == 0) ~"passenger",
    (forms_active_contact == "no" & pharma_important == 0) ~ "bystander",
    TRUE ~ NA_character_
  )) %>%
  mutate(classification = factor(classification, levels = c("bystander", "passenger", "modulator","disconnected driver", "connected driver" )))

# 23 connected drivers, 18 disconnected drivers, 41 modulators, 35 passengers, 278 bystanders
# compare with Figure 5A
# note: if the classification is NA, the mutant is either a low expressor or S74A, which was excluded from classification as efficacy was increased.
residue_classification %>%
  group_by(classification) %>%
  count()


# alternative to forms_active_contacts above, including inactive contacts
# data can be used to generate Figure 4 C/D
b2AR_Gs_data_inactive_active_contacts <- b2AR_Gs_data %>%
  mutate(active_contact = case_when(
    (amino_acid %in% c(contacts_summarised$ResNum1[contacts_summarised$only_in == "3SN6"], contacts_summarised$ResNum2[contacts_summarised$only_in == "3SN6"])) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(inactive_contact = case_when(
    (amino_acid %in% c(contacts_summarised$ResNum1[contacts_summarised$only_in == "2RH1"], contacts_summarised$ResNum2[contacts_summarised$only_in == "2RH1"])) ~ 1,
    TRUE ~ 0
  )) %>%
  right_join(resolved_in_2RH1_3SN6, by = "amino_acid")

# compare Figure 4C
# note that mutation S74A is excluded from the count (listed as NA in "pharma_important")
b2AR_Gs_data_inactive_active_contacts %>%
  filter(expression > 25) %>%
  group_by(pharma_important, active_contact) %>%
  count()

# compare Figure 4D
# note that mutation S74A is excluded from the count (listed as NA in "pharma_important")
b2AR_Gs_data_inactive_active_contacts %>%
  filter(expression > 25) %>%
  group_by(pharma_important, inactive_contact) %>%
  count()


write_csv(residue_classification, "data/residue_classification.csv")
