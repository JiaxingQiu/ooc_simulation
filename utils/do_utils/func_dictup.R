dict_json2df <- function(dict) {
  # return a datafrqme format dictionary with key as rowname dict elememts as colname
  dict_df <- data.frame()
  for (var in names(dict)){
    # name
    var_name <- var
    # label
    var_label <- dict[[var]]$label
    # unique or not
    var_unique <- ifelse(is.null(dict[[var]]$unique_per_sbj), FALSE, dict[[var]]$unique_per_sbj)
    # type
    var_type <- names(dict[[var]])[names(dict[[var]])%in%c("numeric","factor")]
    if (startsWith(var_name,"__")){
      var_type<- "key"
    }
    # unit
    if(var_type=="numeric"){
      var_type<- "num"
      var_unit <- dict[[var]]$numeric$unit
    }else if(var_type=="factor"){
      var_type <- "fct"
      var_unit <- "tag01"
    }else if(var_type=="key"){
      var_unit <- ""
    }else{
      var_unit <- "dictionary error"
      message(paste0("Error prepare dictionary for variable ", var))
    }
    # mlrole
    if (var_type!="key"){
      if ("input" %in% names(dict[[var]]) ){
        if (dict[[var]]$input){
          var_mlrole <- "input"
        }
      }
      if ("output" %in% names(dict[[var]]) ){
        if (dict[[var]]$output){
          var_mlrole <- "output"
        }
      }
    } else {
      var_mlrole <- ""
    }
    dict_df <- bind_rows(dict_df, data.frame(varname=var_name, label=var_label, type=var_type, unit=var_unit, mlrole=var_mlrole, unique_per_sbj=var_unique, stringsAsFactors = FALSE))
  }
  
  # expand factor variables with their levels
  for (var in dict_df$varname[which(dict_df$type=="fct")]){
    fct_df <- data.frame(varname=var, varname_levels=paste0(var,"___",names(dict[[var]]$factor$levels)))
    fct_df <- merge(fct_df, dict_df, all.x=TRUE)
    fct_df <- fct_df[,setdiff(colnames(fct_df), "varname")]
    colnames(fct_df)[which(colnames(fct_df)=="varname_levels")] <- "varname"
    dict_df <- bind_rows(dict_df, fct_df)
  }
  rownames(dict_df) <- dict_df$varname
  
  return(dict_df)
}



get.dict <- function(data, attr_list=c() ){
  if (is.null(attr_list)){
    # get union attr from each column
    attr_list <- names( attributes( data[,1] ) )
    if (length(colnames(data)) >= 2){
      for (i in 2:length(colnames(data))){
        attr_list_new <- names( attributes( data[,i] ) )
        attr_list <- union(attr_list, attr_list_new)
      }
    }
  }
  attr_list <- setdiff(attr_list,c("class","levels"))
  dict_df <- data.frame( matrix( ncol = length(attr_list), nrow = length(colnames(data)) ) )
  colnames(dict_df) <- attr_list
  rownames(dict_df) <- colnames(data)
  for (varname in colnames(data)){
    for (attrname in attr_list){
      dict_df[varname,attrname] <- as.character(ifelse(is.null(attr(data[,varname],attrname)), "", attr(data[,varname],attrname)))
    }
  }
  dict_df$varname <- rownames(dict_df) 
  if(!"label" %in% colnames(dict_df) ){
    dict_df$label <- dict_df$varname
  }
  if(!"unit" %in% colnames(dict_df) ){
    dict_df$unit <- ""
  }
  if(!"type" %in% colnames(dict_df) ){
    dict_df$type <- ""
    for(varname in colnames(data)){
      # automated assign fct and tag01 to certain variables
      if(dplyr::n_distinct(as.character(data[complete.cases(data[,varname]),varname]) )<5){
        dict_df$type[which(dict_df$varname == varname)] <- "fct"
        if(all(unique(as.numeric(as.character(data[complete.cases(data[,varname]),varname]))) %in% c(0,1) )){
          dict_df$unit[which(dict_df$varname == varname)] <- "tag01"
        }
      }else{
        dict_df$type[which(dict_df$varname == varname)] <- "num"
      }
    }
  }
  if("unit" %in% colnames(dict_df)){
    if(!"unit_label" %in% colnames(dict_df) ){
      dict_df$unit_label <- ""
    }
    dict_df$unit_label[which(dict_df$unit=="tag01")] <- "1=Yes; 0=No" 
  }
  if(!"unique_per_sbj" %in% colnames(dict_df) ){
    dict_df$unique_per_sbj <- "FALSE"
  }
  if(!"source_file" %in% colnames(dict_df) ){
    dict_df$source_file <- paste0("Created at ",Sys.time())
  }
  return(dict_df)
}


remove.dict <- function(data){
  for (col in colnames(data)) attributes(data[[col]]) <- NULL
  return(data)
}

assign.dict <- function(data, dict_new, multi_assign=FALSE, overwrite=TRUE){
  # ---- inputs ----
  # data: data frame object
  # dictionary data frame object
  # overwrite: whether or not overwrite attributes that already exist in the dataframe (default true)
  # multi_assign: whether or not rownames in dict_df match with multiple columns in the dataframe
  
  # ---- output ----
  # dataframe object with updated attributes
  
  #  ---- make valid dict_new
  stopifnot(!is.null(dict_new$varname))
  dict_new$varname[which(dict_new$varname=="")] <-NA
  dict_new <- dict_new %>% filter(!is.na(varname)) %>% as.data.frame()
  dict_new <- dplyr::distinct(dict_new)
    
  # colnames(data) <- gsub("[^[:alnum:]]+","_",colnames(data)) # save for later
  # dict_new$varname <- gsub("[^[:alnum:]]+","_",dict_new$varname)
  
  colnames(data) <- gsub("[^[:alnum:]]","_",colnames(data))
  dict_new$varname <- gsub("[^[:alnum:]]","_",dict_new$varname)
  
  if ("varname_dict" %in% colnames(dict_new)) dict_new$varname_dict <- gsub("[^[:alnum:]]","_",dict_new$varname_dict)
  
  
  if(multi_assign) {
    dict_new_key <- data.frame(stringsAsFactors = FALSE)
    colnames(dict_new)[which(colnames(dict_new)=="varname")] <- "varname_dict"
    for (var in dict_new$varname_dict) {
      for(col in colnames(data)){
        if(startsWith(col,var)|endsWith(col,var)|startsWith(col,paste0("X",var))){
          dict_new_key <- bind_rows(dict_new_key,data.frame(varname_dict=var, varname=col , stringsAsFactors = FALSE))
        }
      }
    }
    dict_new <- merge(dict_new, dict_new_key)
  }
  
  rownames(dict_new) <- as.character( dict_new$varname )
  
  # keep old dictionary
  dict_old <- get.dict(data)
  
  # clean missing value in new or old dictionary
  for (col in colnames(dict_new)){
    dict_new[,col] <- as.character(dict_new[,col])
    dict_new[which( is.na(dict_new[,col])|dict_new[,col]=="NA" ),col] <- ""
  }
  for (col in colnames(dict_old)){
    dict_old[,col] <- as.character(dict_old[,col])
    dict_old[which( is.na(dict_old[,col])|dict_old[,col]=="NA" ),col] <- ""
  }
  
  dict_df <- data.frame( matrix(data="", ncol = length(union(colnames(dict_old),colnames(dict_new))), nrow = length(union(rownames(dict_old),rownames(dict_new))) ) , stringsAsFactors = FALSE)
  colnames(dict_df) <- union(colnames(dict_old),colnames(dict_new))
  rownames(dict_df) <- union(rownames(dict_old),rownames(dict_new))
  
  if (!overwrite){
    for (var in rownames(dict_new)){
      for (attr in colnames(dict_new)){
        if(dict_new[var, attr]!=""){
          dict_df[var, attr] <- dict_new[var, attr]
        }
      }
    }
    for (var in rownames(dict_old)){
      for (attr in colnames(dict_old)){
        if(dict_old[var, attr]!=""){
          dict_df[var, attr] <- dict_old[var, attr]
        }
      }
    }
  }else{
    for (var in rownames(dict_old)){
      for (attr in colnames(dict_old)){
        if(dict_old[var, attr]!=""){
          dict_df[var, attr] <- dict_old[var, attr]
        }
      }
    }
    for (var in rownames(dict_new)){
      for (attr in colnames(dict_new)){
        if(dict_new[var, attr]!=""){ 
          dict_df[var, attr] <- dict_new[var, attr]
        }
      }
    }
  }
  
  # find vars to assign attributes
  vars2assign <- intersect(colnames(data), rownames(dict_df))
  for (varname in vars2assign){
    for (attrname in colnames(dict_df) ){
      attr( data[,varname], attrname ) <- as.character( dict_df[varname, attrname] )
    }
  }
  
  return(data)
}



merge_with_dict <- function(x,y, all=FALSE, all.x=all, all.y=all){
  dict.x <- get.dict(x)
  dict.y <- get.dict(y)
  df <- merge(x, y, all.x=all.x,all.y=all.y,all=all)
  if (dim(dict.y)[1]*dim(dict.y)[2]!=0){
    df <- assign.dict(df, dict.y)
  }
  if (dim(dict.x)[1]*dim(dict.x)[2]!=0){
    df <- assign.dict(df, dict.x)
  }
  return(df)
}



dict.tag2fct <- function(data, revlist=NULL, use_label=FALSE, fillna="None_level"){
  
  dict_df <- get.dict(data)
  stopifnot("varname_dict"%in%colnames(dict_df))
  #stopifnot("mlrole"%in% colnames(dict_df))
  stopifnot("tag01"%in%unique(dict_df$unit))
  data <- assign.dict(data, dict_df)
  
  tag_dict <- dict_df[which(dict_df$type=="fct" & dict_df$unit=="tag01" ),]# & dict_df$mlrole%in%c("input","output","cluster")
  revlist_dict <- unique(tag_dict$varname_dict)
  if(is.null(revlist)){
    revlist <- revlist_dict
  }else{
    revlist <- intersect(revlist, revlist_dict) 
  }
  
  # create reversed factor columns
  for (fct_col in revlist){
    data[,paste0(fct_col,"___rev")] <- fillna # defualt new column
    for (tag_col in tag_dict$varname[which(tag_dict$varname_dict==fct_col)]) {
      print(paste0("--- reverse tag_col ", tag_col, " to fct_col ", paste0(fct_col,"___rev"), " ---"))
      if (any(unique(data[,tag_col])==1)){
        data[which(data[,tag_col]==1),paste0(fct_col,"___rev")] <- ifelse(use_label,  paste0(gsub("[^[:alnum:]]","_",attr(data[,tag_col], "label")),"_", gsub(fct_col,"",tag_col)) , gsub(fct_col,"",tag_col) )
        
      }else if (any(unique(data[,tag_col])%in%c("Yes","yes")  )){
        data[which(data[,tag_col]%in%c("Yes","yes") ),paste0(fct_col,"___rev")] <- ifelse(use_label,paste0(gsub("[^[:alnum:]]","_",attr(data[,tag_col], "label")),"_", gsub(fct_col,"",tag_col)) ,  gsub(fct_col,"",tag_col))
        
      } else {
        message(paste0(" --- unrecoganized level in tag_col ---", tag_col))
      }
    }
    
    for(att in names(attributes(data[,tag_col])) ){
      attr(data[,paste0(fct_col,"___rev")], att) <- attr(data[,tag_col], att)
    }
    attr(data[,paste0(fct_col,"___rev")], "varname_dict") <- fct_col
    attr(data[,paste0(fct_col,"___rev")], "varname") <- paste0(fct_col,"___rev")
    attr(data[,paste0(fct_col,"___rev")], "label") <- fct_col
    attr(data[,paste0(fct_col,"___rev")], "type") <- "fct"
    attr(data[,paste0(fct_col,"___rev")], "unit") <- ""
  }
  
  dict_df_new = get.dict(data)
  dictup_data <- assign.dict(data, dict_df_new)
  
  return(dictup_data)
  
}






# assign.dict <- function(data, dict){
#   # data: data frame object
#   # dictionary data frame object
#   for (var in as.character(dict$name)){
#     if( var %in% colnames(data)){
#       attr(data[,var],"label") <- ifelse(is.na(dict$label[which(dict$name==var)]),"",as.character(dict$label[which(dict$name==var)]))
#       attr(data[,var],"type") <- ifelse(is.na(dict$type[which(dict$name==var)]),"",as.character(dict$type[which(dict$name==var)]))
#       attr(data[,var],"unit") <- ifelse(is.na(dict$unit[which(dict$name==var)]),"",as.character(dict$unit[which(dict$name==var)]))
#     }
#     if( paste0(var,".factor") %in% colnames(data)){
#       attr(data[,paste0(var,".factor")],"label") <- ifelse(is.na(dict$label[which(dict$name==var)]),"",as.character(dict$label[which(dict$name==var)]))
#       attr(data[,paste0(var,".factor")],"type") <- ifelse(is.na(dict$type[which(dict$name==var)]),"",as.character(dict$type[which(dict$name==var)]))
#       attr(data[,paste0(var,".factor")],"unit") <- ifelse(is.na(dict$unit[which(dict$name==var)]),"",as.character(dict$unit[which(dict$name==var)]))
#     }
#   }
#   for (var in colnames(data)[which(grepl("___tag",colnames(data)))]){
#     attr(data[,var],"type") <- "fct"
#     attr(data[,var],"unit") <- "tag01"
#   }
#   return(data)
# }



# get.dict <- function(data){
#   # input:
#   #    data: data.frame object with assigned attributes
#   #    dict_path: file path to save dictionary csv 
#   dict <- data.frame(name=c(),label=c(),type=c(),unit=c())
#   for (var in colnames(data)){
#     var_name <- var
#     var_label <- as.character(ifelse(is.null(attr(data[,var],"label")), "", attr(data[,var],"label")))
#     var_type <- as.character(ifelse(is.null(attr(data[,var],"type")), "", attr(data[,var],"type")))
#     var_unit <- as.character(ifelse(is.null(attr(data[,var],"unit")), "", attr(data[,var],"unit")))
#     dict <- bind_rows(dict,data.frame(name=var_name, label=var_label, type=var_type, unit=var_unit))
#   }
#   return(dict)
# }

