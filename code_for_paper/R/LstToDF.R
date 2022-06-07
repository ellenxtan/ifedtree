## Transform nested lists to dataframe (data.table)
LstToDF <- function(lst, df_name) {
    out_lst <- purrr::map(lst, df_name)
    out_df <- purrr::map_dfr(out_lst, rbind)  # row-bind
    return(out_df)
}
