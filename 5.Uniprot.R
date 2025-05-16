#Load necessary libraries
library(tidyverse)

#Read and process the tsv files
uniprot <- read_tsv('idmapping_2024_11_13.tsv', col_types = cols(.default = "c"))
write.csv(uniprot,"uniprot.csv")

#Allergens
allergens <- uniprot %>%
  filter(!is.na(`Allergenic Properties`)) %>%
  select(`Entry Name`, `Protein names`, `Allergenic Properties`, `Organism`) %>%
  group_by(`Allergenic Properties`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(allergens,"Allergens.csv")

#Disease
disease <- uniprot %>%
  filter(!is.na(`Involvement in disease`)) %>%
  select(`Entry Name`, `Protein names`, `Involvement in disease`, `Organism`) %>%
  group_by(`Involvement in disease`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(disease,"Disease.csv")

#Virus hosts
virushosts <- uniprot %>%
  filter(!is.na(`Virus hosts`)) %>%
  select(`Entry Name`, `Protein names`, `Virus hosts`, `Organism`) %>%
  group_by(`Virus hosts`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(virushosts,"Virushosts.csv")

#Post-translational modification
ptm <- uniprot %>%
  filter(!is.na(`Post-translational modification`)) %>%
  select(`Entry Name`, `Protein names`, `Post-translational modification`, `Organism`) %>%
  group_by(`Post-translational modification`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(ptm,"ptm.csv")

#Kinetics
kinetics <- uniprot %>%
  filter(!is.na(`Kinetics`)) %>%
  select(`Entry Name`, `Protein names`, `Kinetics`, `Organism`) %>%
  group_by(`Kinetics`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(kinetics,"Kinetics.csv")

#Catalytic activity
Catalytic_activity <- uniprot %>%
  filter(!is.na(`Catalytic activity`)) %>%
  select(`Entry Name`, `Protein names`, `Catalytic activity`, `Organism`) %>%
  group_by(`Catalytic activity`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Catalytic_activity,"Catalytic_activity.csv")

#Zinc finger
zf <- uniprot %>%
  filter(!is.na(`Zinc finger`)) %>%
  select(`Entry Name`, `Protein names`, `Zinc finger`, `Organism`) %>%
  group_by(`Zinc finger`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(zf,"zf.csv")

#Cofactor
Cofactor <- uniprot %>%
  filter(!is.na(`Cofactor`)) %>%
  select(`Entry Name`, `Protein names`, `Cofactor`, `Organism`) %>%
  group_by(`Cofactor`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Cofactor,"Cofactor.csv")

#DNA binding
DNAbinding <- uniprot %>%
  filter(!is.na(`DNA binding`)) %>%
  select(`Entry Name`, `Protein names`, `DNA binding`, `Organism`) %>%
  group_by(`DNA binding`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(DNAbinding,"DNA_binding.csv")

#Toxic dose
Toxicdose <- uniprot %>%
  filter(!is.na(`Toxic dose`)) %>%
  select(`Entry Name`, `Protein names`, `Toxic dose`, `Organism`) %>%
  group_by(`Toxic dose`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Toxicdose,"Toxic_dose.csv")

#Pharmaceutical use
Pharmaceutical <- uniprot %>%
  filter(!is.na(`Pharmaceutical use`)) %>%
  select(`Entry Name`, `Protein names`, `Pharmaceutical use`, `Organism`) %>%
  group_by(`Pharmaceutical use`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Pharmaceutical,"Pharmaceutical.csv")

#Transmembrane
Transmembrane <- uniprot %>%
  filter(!is.na(`Transmembrane`)) %>%
  select(`Entry Name`, `Protein names`, `Transmembrane`, `Organism`) %>%
  group_by(`Transmembrane`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Transmembrane,"Transmembrane.csv")

#Intramembrane
Intramembrane <- uniprot %>%
  filter(!is.na(`Intramembrane`)) %>%
  select(`Entry Name`, `Protein names`, `Intramembrane`, `Organism`) %>%
  group_by(`Intramembrane`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Intramembrane,"Intramembrane.csv")

#Signal peptide
Signalpeptide <- uniprot %>%
  filter(!is.na(`Signal peptide`)) %>%
  select(`Entry Name`, `Protein names`, `Signal peptide`, `Organism`) %>%
  group_by(`Signal peptide`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Signalpeptide,"Signal_peptide.csv")

#Motif
Motif <- uniprot %>%
  filter(!is.na(`Motif`)) %>%
  select(`Entry Name`, `Protein names`, `Motif`, `Organism`) %>%
  group_by(`Motif`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Motif,"Motif.csv")

#Protein families
Proteinfamilies <- uniprot %>%
  filter(!is.na(`Protein families`)) %>%
  select(`Entry Name`, `Protein names`, `Protein families`, `Organism`) %>%
  group_by(`Protein families`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Proteinfamilies,"Protein_families.csv")

#DrugBank
DrugBank <- uniprot %>%
  filter(!is.na(`DrugBank`)) %>%
  select(`Entry Name`, `Protein names`, `DrugBank`, `Organism`) %>%
  group_by(`DrugBank`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(DrugBank,"DrugBank.csv")

#Allergome
Allergome <- uniprot %>%
  filter(!is.na(`Allergome`)) %>%
  select(`Entry Name`, `Protein names`, `Allergome`, `Organism`) %>%
  group_by(`Allergome`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Allergome,"Allergome.csv")

#Antibodypedia
Antibodypedia <- uniprot %>%
  filter(!is.na(`Antibodypedia`)) %>%
  select(`Entry Name`, `Protein names`, `Antibodypedia`, `Organism`) %>%
  group_by(`Antibodypedia`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(Antibodypedia,"Antibodypedia.csv")

#GeneTree
GeneTree <- uniprot %>%
  filter(!is.na(`GeneTree`)) %>%
  select(`Entry Name`, `Protein names`, `GeneTree`, `Organism`) %>%
  group_by(`GeneTree`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(GeneTree,"GeneTree.csv")

#eggNOG
eggNOG <- uniprot %>%
  filter(!is.na(`eggNOG`)) %>%
  select(`Entry Name`, `Protein names`, `eggNOG`, `Organism`) %>%
  group_by(`eggNOG`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(eggNOG,"eggNOG.csv")

#PlantReactome
PlantReactome <- uniprot %>%
  filter(!is.na(`PlantReactome`)) %>%
  select(`Entry Name`, `Protein names`, `PlantReactome`, `Organism`) %>%
  group_by(`PlantReactome`) %>%
  summarise_all(list) %>%
  ungroup()%>%
  unnest(cols = everything())
write.csv(PlantReactome,"PlantReactome.csv")

#Merge DrugBank and Disease dataframes on the 'Entry Name' column
drug_disease_interaction <- merge(x = DrugBank, y = disease, by = "Protein names", suffixes = c("_Drug", "_Disease"))%>%
  select(`Protein names`, DrugBank, `Involvement in disease`)
write.csv(drug_disease_interaction,"drug_disease_interaction.csv")

#counts for the drugs dataset
drug_counts <- drug_disease_interaction %>%
  separate_rows(DrugBank, sep = ";") %>%
  mutate(DrugBank = str_trim(DrugBank),
         DrugBank = str_remove(DrugBank, "\\.$")) %>%
  filter(!str_starts(DrugBank, "DB")) %>%
  count(DrugBank, sort = TRUE) %>% 
  filter(n>5)
write.csv(drug_counts,"drug_counts.csv")


