# DeterminantsSignalingGPCR

Repository for the manuscript "Molecular determinants of ligand efficacy and potency in GPCR signaling".
Contains the relevant code to understand and assess the conclusions of the manuscript.

Paths might need to be changed in scripts to match file location and requirements of the operating system.
Scripts import data made available in the data folder of this repository. All scripts can be run using R or RStudio.

## Data annotation

The data_annotation.R script imports raw and normalised Gs signalling data for 412 A/G mutations of the beta2-adrenegic receptor in response to adrenaline. For more information on the methods and experiment design, please check the manuscript.

The script further imports and cleans Evolutionary Trace (ET) Scores obtained from http://evolution.lichtargelab.org/gpcr. This step requires GPCR residue-number tables that are provided in the data folder.

Structures for PDB ID 2RH1 and 3SN6 from https://www.rcsb.org/ were aligned and .pdb files were exported for further analysis in this script. We have extracted x, y, and z coordinates and dihedral angles for both structures. Accessible surface area data were calculated using dssp https://swift.cmbi.umcn.nl/gv/dssp/ and are provided in the data folder.

Code for combination and further annotation of the data are provided in this script. The script illustrates how Table S1 of the manuscript was assembled. The final table is exported as a .csv file and is required for the subsequent scripts.

## Network

The network.R script illustrates the processing of atom-atom contacts for PDB IDs 2RH1 and 3SN6 obtained from https://pca.mbgroup.bio/index.html. First, the contact tables as downloaded from the website are imported and cleaned. We then list the residues resolved in both structures, combine the contact tables and show how we reached the numbers quoted in Figure 1C of the manuscript. We then list all residues where mutation led to decreased potency and/or efficacy as well as their contacts and illustrate the analysis steps to reach the network (cf. Figure 5).

We illustrate the residue classification, which is exported as a .csv file into the data folder. This file is required for the analysis of the allosteric modulators as shown in Figure 6.

In addition, we provide the counts of pharmacologically important residues that do or don not form an inactive or active state-specific contact as listed and plotted in Figure 4 C/D.

## Ligand binding

The ligand_binding.R script imports contact tables for the orthosteric ligand adrenaline and the allosteric modulators 6FA (positive allosteric modulator, PAM) and AS408 (negative allosteric modulator, NAM) and illustrates the processing of the data. Results can be compared with the data shown in Figure 3B (adrenaline), Figure 6E (NAM) and Figure 6F (PAM).

## G protein binding

The Gprotein_binding.R script imports contacts between the G protein (alpha and beta subunits) and the receptor, using the active, G protein-bound structure 3SN6. Furthermore, Common G protein numbering (CGN) is imported from the CGN server, https://www.mrc-lmb.cam.ac.uk/CGN/. Finally, we show the data underlying Figure 3E.
