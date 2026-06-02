########
#testing
########

#libraries 
source("/datos/home/lorepaty/scripts/functions_mi.R")


#paths

path_micro <- "/datos/ihmp/2026-05-29-normalizado_filtrado/non_micro_filtered_CLR_normalized.csv"
path_tx <- "/datos/ihmp/2026-05-29-normalizado_filtrado/non_tx_CPM_log_normalized.csv"
path_output <- "/datos/ihmp/2026-05-29-normalizado_filtrado/results_2026_05_29/non_MI_norm.csv"

#read data 

micro3 <- as.data.frame(readr::read_csv(path_micro))
tx3 <- as.data.frame(readr::read_csv(path_tx))




#micro3 <- micro[1:10, 1:10]
#tx3        <- tx[1:10, 1:10]
#used for testing


#discretizing 

tempus <- proc.time()
d.micro <- par_discretizer(micro3, korez = 6)
tempus <- proc.time() - tempus
print(tempus)

tempus <- proc.time()
d.tx <- par_discretizer(tx3, korez = 6)
tempus <- proc.time() - tempus
print(tempus)

#mi calculating 

tempus <- proc.time()
mirXrna <- par_mi_calc(sources = d.micro, 
                       targets = d.tx, 
                       korez = 6)
tempus <- proc.time() - tempus
print(tempus)

mi_matrix <- bind_rows(!!!mirXrna, #explicit splicing
                       .id = "micro/tx")

write_csv(mi_matrix, path_output)