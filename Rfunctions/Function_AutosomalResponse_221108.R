AutosomalResponse <- function(dataset){
  subset(dataset, dataset$chr != "chrX" & dataset$chr != "chrY")
}

