#' Theory behind association studies

#' PLINK files
#' PED file
ped <- data.frame(FamilyID = "FAM001",
                  IndividualID = 1:7,
                  PaternalID = c(0,0,0,1,1,1,5),
                  MaternalID = c(0,0,0,2,2,3,6),
                  Sex = c(1,2,2,1,1,2,2),
                  Phenotype = c(0,1,1,1,1,0,0),
                  snp1a = c("A","B","A","A","A","A","A"),
                  snp1b = c("A","B","B","B","B","A","A")
                  )
write.table(ped, "data/assoc_01.ped",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
map <- data.frame(chromosome = 1,
                  snp = "snp1",
                  cm = 0,
                  bp = 1000
                  )
write.table(map, "data/assoc_01.map",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
