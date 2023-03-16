df_to_workbook <- function(file_name, data, have_rownames=TRUE) {
  # save multiple data frame to a single excel spreadsheet
  #@ file_name: string, the name of the excel file
  #@ data: a named vector or list. name is the worksheet name,
  #        the value is the data frame
  
  # Create Workbook object
  wb <- createWorkbook()
  
  # Add each dataframe as a worksheet using map2()
  map2(names(data), data, function(name, df) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, df, 
              startCol = 1, startRow = 1, rowNames = have_rownames)
  })
  
  # Save the workbook to the specified file name
  saveWorkbook(wb, file_name)
}

get_res_from_dds <- function(dds, contrast, res_filter=NULL ) {
  # get the DEG result based on contrast matrix and return the 
  # filterd result if the filter is provided
  
  if (!inherits(dds, "DESeqDataSet")){
    stop("The object must be a DESeqDataSet!")
  }
  
  
  res = results(dds, contrast = contrast)
  res = as.data.frame(res)
  
  if(is.null(res_filter)){
    return (res)
  } else{
    return (subset(res, 
                   subset = eval(parse(text = res_filter)))
    )
  } 
  
}
extract_immune_gene_lst <- function(gene_name_lst){
  # get the inflammation and immune response gene list
  pattern =  "^IL|^CCR|^CXCR|^PTGS|^NLRP3|^CD|^IFN|^TLR|^NF"
  immune_genes = grep(pattern = pattern, gene_name_lst, value = T)
  
  return (immune_genes)
}
