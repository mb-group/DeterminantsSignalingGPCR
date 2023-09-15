# TO RUN: download the data folder and change the location in the setwd() command to the parent directory of the data folder
# Windows users will need to change "data/" to "data\\" throughout the script

# Inputs: the script requires basic signalling data, cut-offs, Evolutionary Trace Scores, tables with residue and GPCRdb numbers
# Evolutionary Trace Scores were obtained from: http://evolution.lichtargelab.org/gpcr
# Here, we show how the different inputs are cleaned/processed 
# We use bio3d functions for structure analysis: http://thegrantlab.org/bio3d/
# Output: a single table containing most of the measures listed in Table S1. The table is required for running additional scripts. 


# install required packages
packages <- c("tidyverse", "bio3d")
install.packages(setdiff(packages, rownames(installed.packages())))

# load and attach add-on packages
library(tidyverse)
library(bio3d)

# change the working directory as needed. The working directory should be the parent folder of the 'data' folder. 
# please note that paths require double backslashes when running scripts on Windows.

setwd("/Users/xyz")

# import signalling data & cut-offs
b2AR_Gs_cutoffs <- read_csv("data/b2AR_Gs_cutoffs.csv")
b2AR_Gs_raw <- read_csv("data/b2AR_Gs_raw.csv")

# Evolutionary Trace
# Data used for Figures 6, S13

# import Evolutionary Trace data
ET_adrenoceptor_import <- read_table("data/Adrenoceptor.txt", col_names = TRUE)
ET_aminergic_import <-  read_table("data/All_Amine.txt", col_names = TRUE)
ET_classA_import <- read_table("data/ClassA_5105.txt", col_names = TRUE)

# import residue table for Rhodopsin (amino acid numbers in ET data files are Rhodopsin numbers)
# tables with GPCRdb numbers were downloaded from gpcrdb.org, gaps were removed, and positions without GPCRdb number were annotated
residue_table_rho <- read_csv("data/bovine_rhodopsin_residue_table.csv")
residue_table_hb2AR <- read_csv("data/hb2AR_residue_table_clean.csv")

# function for cleaning ET import files
clean_ET <- as_mapper(~.x %>%
                        mutate_at(vars(alignment, residue, coverage, variability_n, Score), as.numeric) %>%
                        left_join(residue_table_rho %>% select(amino_acid, GPCRdb), by = c("residue" = "amino_acid")) %>%
                        # ECL2 residues cannot be matched to b2AR correctly
                        filter(GPCRdb != "ECL2"& GPCRdb != "ECL3") %>%
                        left_join(residue_table_hb2AR %>% select(b2AR_amino_acid = residue, GPCRdb), by = "GPCRdb"))

# clean up imported ET tables
ET_adrenoceptor <- clean_ET(ET_adrenoceptor_import)
ET_aminergic <- clean_ET(ET_aminergic_import)
ET_classA <- clean_ET(ET_classA_import)



# Structures


# coordinates
# 2RH1 and 3SN6 were aligned in pymol and the new coordinates exported
# coordinates used for Figures 5, S4, S13
# distances used for Figures 4, S1

b2AR_3SN6_new <- read.pdb("data/3SN6_new_coordinates.pdb")
b2AR_3SN6_new_coordinates <- b2AR_3SN6_new$atom %>%
  filter(elety == "CA")%>%
  select(resno, x, y, z)
b2AR_3SN6_new_coordinates_long <- pivot_longer(b2AR_3SN6_new_coordinates, 2:4, names_to = "coord")

b2AR_2RH1_new <- read.pdb("data/2RH1_new_coordinates.pdb")
b2AR_2RH1_new_coordinates <- b2AR_2RH1_new$atom %>%
  filter(elety == "CA")%>%
  select(resno, x, y, z)
b2AR_2RH1_new_coordinates_long <- pivot_longer(b2AR_2RH1_new_coordinates, 2:4, names_to = "coord")

b2AR_3SN6_2RH1_coordinates <- inner_join(b2AR_3SN6_new_coordinates, b2AR_2RH1_new_coordinates, suffix = c("_3SN6", "_2RH1"), by = "resno")%>%
  mutate(distance = sqrt((x_3SN6-x_2RH1)^2+(y_3SN6-y_2RH1)^2+(z_3SN6-z_2RH1)^2))

# torsion angles
# Data used for Figure S7

torsion_3SN6 <- torsion.pdb(b2AR_3SN6_new)
torsion_2RH1 <- torsion.pdb(b2AR_2RH1_new)

angles_3SN6 <- torsion_3SN6$tbl %>%
  as_tibble() %>%
  mutate(rowname = rownames(torsion_3SN6$tbl), 
         resno = as.numeric(str_extract(rownames(torsion_3SN6$tbl), "\\d{1,3}")), 
         omega = torsion_3SN6$omega) %>%
  filter(rowname != "1601.R.P0G") %>%
  select(rowname, resno, everything())
angles_2RH1 <- torsion_2RH1$tbl %>% 
  as_tibble() %>% 
  mutate(rowname = rownames(torsion_2RH1$tbl), 
         resno = as.numeric(str_extract(rownames(torsion_2RH1$tbl), "\\d{1,3}")), 
         omega = torsion_2RH1$omega) %>%
  filter(resno < 399) %>%
  select(rowname, resno, everything())

angles <- angles_3SN6 %>%
  inner_join(angles_2RH1, by = "resno", suffix = c("_3SN6", "_2RH1")) %>%
  mutate(delta_phi = round((phi_3SN6-phi_2RH1+180)%%360-180, 1), 
         delta_psi = round((psi_3SN6-psi_2RH1+180)%%360-180, 1), 
         delta_omega = round((omega_3SN6-omega_2RH1+180)%%360-180, 1))

# import accessible surface area data. Data were computed using dssp https://swift.cmbi.umcn.nl/gv/dssp/. 
# Since the calculations requires downloading dssp, we only provide the cleaned output here. dssp was run on 2RH1 and 3SN6.
# Data used for Figure 2C.

dssp_ASA <- read_csv("data/ASA_dssp_2RH1_3SN6.csv")


# combining data

b2AR_Gs_data <- b2AR_Gs_raw %>%
  mutate(efficacy = case_when(
    expression < 25 ~ NA_character_,
    amplitude < b2AR_Gs_cutoffs$min_amp ~ "reduced",
    amplitude > b2AR_Gs_cutoffs$max_amp ~ "increased",
    amplitude >= b2AR_Gs_cutoffs$min_amp ~ "WT-like",
    TRUE ~ "error"
  )) %>%
  mutate(potency = case_when(
    expression < 25 ~ NA_character_,
    curation == 2 ~ "reduced",
    EC50 > b2AR_Gs_cutoffs$max_logEC50 ~ "reduced",
    EC50 < b2AR_Gs_cutoffs$min_logEC50 ~ "increased",
    EC50 <= b2AR_Gs_cutoffs$max_logEC50 ~ "WT-like",
    TRUE ~ "error"
  )) %>%
  mutate(pharma_important = case_when(
    expression < 25 ~ NA_real_,
    (potency == "reduced" | efficacy == "reduced") ~ 1,
    (potency == "WT-like" & efficacy == "WT-like") ~ 0,
    efficacy == "increased" ~ NA_real_
  )) %>%
  mutate(x50 = if_else(str_detect(GPCRdb, "^[1-8]x50"), 1, 0)) %>%
  mutate(motif = case_when(
    GPCRdb %in% c("3x49", "3x50", "3x51") ~ "DRY",
    GPCRdb %in% c("6x47", "6x48", "6x50") ~ "CWxP",
    GPCRdb %in% c("7x49", "7x50", "7x53") ~ "NPxxY",
    GPCRdb %in% c("3x40", "5x50", "6x44") ~ "PIF",
    TRUE ~ NA_character_
  )) %>%
  left_join(ET_adrenoceptor %>% select(b2AR_amino_acid, ET_Score_adrenergic = Score), by = c("amino_acid" = "b2AR_amino_acid")) %>%
  left_join(ET_aminergic %>% select(b2AR_amino_acid, ET_Score_aminergic = Score), by = c("amino_acid" = "b2AR_amino_acid")) %>%
  left_join(ET_classA %>% select(b2AR_amino_acid, ET_Score_classA = Score), by = c("amino_acid" = "b2AR_amino_acid")) %>%
  left_join(b2AR_3SN6_2RH1_coordinates, by = c("amino_acid" = "resno")) %>%
  left_join(angles %>% select(contains("delta"), resno), by = c("amino_acid" = "resno")) %>%
  left_join(dssp_ASA, by = c("amino_acid" = "resno")) %>%
  mutate(delta_acc = acc_2RH1-acc_3SN6)

# 82 pharmacologically important residues
# compare with Figure 5C
# compare main text, subheading 'Conserved motifs only explain a fraction of residues important for signaling'
# compare main text, subheading 'One-fifth of the positions in β2AR are important for adrenaline efficacy, potency, or both'
b2AR_Gs_data %>% 
  filter(pharma_important == 1)

# two-thirds of the pharmacologically important residues (55 of 82 residues) do not map to known functional sites or motifs
# compare main text, subheading 'Conserved motifs only explain a fraction of residues important for signaling'
b2AR_Gs_data %>% 
  filter(pharma_important == 1) %>%
  filter(is.na(motif)) %>%
  # 74 residues left after motif filter
  filter(Gprotein == 0) %>%
  # 65 residues left after Gprotein filter
  filter(adrenaline == 0)
  # 55 residues left after adrenaline binding site filter
  # 55/82 residues => 67%. 


write_csv(b2AR_Gs_data, "data/b2AR_Gs_data.csv")