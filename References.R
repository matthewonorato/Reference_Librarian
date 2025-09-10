#' @title Get paths to annotation files or BED files
get_paths <- function(refs=NULL, file_id=NULL, ref_path=NULL){

    stopifnot("refs cannot be NULL" = !is.null(refs))
    stopifnot("file_id cannot be NULL" = !is.null(file_id))
    stopifnot("ref_path cannot be NULL" = !is.null(ref_path))
    
    file_paths = list.files(path=ref_path,
                            pattern=paste0(".*",file_id,".*"),
                            all.files=FALSE, full.names=TRUE, recursive=TRUE, ignore.case=FALSE)
    
    # Ensure same order for ref names and file paths
    file_pattern = tolower(gsub(" ", ".*", refs, fixed=TRUE))
    file_paths = map_chr(file_pattern, ~ { file_paths[str_detect(tolower(file_paths), .x)] })
    return(file_paths)
}

#' @title Concatenate annotation files or BED files
concat_files <- function(file_paths, ref_names, ref_colors, to_chr=NULL) {
    
    stopifnot(length(file_paths) == length(ref_names),
              length(ref_names) == length(ref_colors))

    bed_rename = c(V1="chr", V2="start", V3="end", V4="feat_id", V5="score", V6="strand")
    
    df = pmap(.l = list(file_path = file_paths, ref_name = ref_names, ref_color = ref_colors),
              .progress = FALSE,
              .f <- function(file_path, ref_name, ref_color) {
                  fread(file_path, data.table=FALSE) %>%  # Read in as dataframe
                  rename_with(~ bed_rename[.x], .cols=any_of(names(bed_rename))) %>%
                  clean_names() %>%
                  mutate(ref_name = ref_name,
                         color = ref_color,
                         seq_id = paste0("[", ref_name, "] ", chr, ":", start, "-", end),
                         length = case_when(
                             grepl("\\.bed", file_path, ignore.case=TRUE) ~ end-start,  # BED is 0-based
                             .default=end-start+1)  # Other files are 1-based, so add 1
                        )}) %>% 
         bind_rows()  # Bind lists and convert to df
  
    if (!is.null(to_chr)) {  # Filter to primaries only
        df = df %>% filter(chr %in% to_chr)
        message("Filtered to chromosomes: ", toString(to_chr)) }
    return(df)
}

#' @title Get unique features
#'
#' @note df must contain columns: feat_id, ref_name
get_unique_feats <- function(df, feat_type=NULL){
    
    feat_type = arg_match(feat_type, c("gene_name", "gene_id", "transcript_id"), error_arg="feat_type")
    feat_type_sym = sym(feat_type)
    
    df_unique = df %>% mutate(feat_id = {{feat_type_sym}}) %>%
        mutate(across(all_of("feat_id"), ~ str_remove(.x, "\\..*"))) %>%  # Remove version numbers
        group_by(ref_name) %>%
        distinct(feat_id, .keep_all=TRUE) %>%  # Remove duplicate features
        filter(!is.na(feat_id) & feat_id != "") %>% ungroup()  # Remove missing features
    return(df_unique)
}

#' @title Get total counts
#'
#' @note df must contain columns: feat_id, ref_name, color
get_totals <- function(df){
    df_total = df %>% group_by(ref_name, color) %>%
        summarize(Count=n_distinct(feat_id), .groups="keep") %>% ungroup() %>%
        mutate(chr="All\nPrimaries")
    return(df_total)
}

#' @title Get total counts per chromosome
#'
#' @note df must contain columns: feat_id, ref_name, color, chr
get_totals_per_chr <- function(df){
    df_total_per_chr = df %>% group_by(ref_name, color, chr) %>%
        summarize(Count=n_distinct(feat_id), .groups="keep") %>% ungroup() %>%
        mutate(ref_name = factor(ref_name, levels=unique(ref_name)),
               chr = factor(chr, levels=primary_chr))
    return(df_total_per_chr)
}

#' @title Get feature types
#'
#' @note df must contain columns: gene_biotype or transcript_biotype, ref_name, color
get_totals_per_type <- function(df, feat_type=NULL, types_to_plot=c("protein_coding")){
    
    feat_type = arg_match(feat_type, c("gene_name", "gene_id", "transcript_id"), error_arg="feat_type")

    if (feat_type == "transcript_id") {
        df = df %>% mutate(feat_biotype = transcript_biotype)
    } else {
        df = df %>% mutate(feat_biotype = gene_biotype) }
    
    types_to_plot = arg_match(types_to_plot, unique(df$feat_biotype), error_arg="types_to_plot")
    
    df_types = df %>% filter(feat_biotype %in% types_to_plot) %>%
        group_by(ref_name, color, feat_biotype) %>%
        summarize(Count=n(), .groups="keep") %>% ungroup() %>%
        mutate(chr="All\nPrimaries", feat_biotype = factor(feat_biotype, levels=types_to_plot))
    return(df_types)
}

#' @title Create dumbbell plot
#'
#' @description Create dumbbell plot of feature counts across any number of references
plot_ref_counts <- function(df, title_text=NULL, subtitle_text=NULL, out_file=NULL, w=6, h=3) {

    stopifnot("out_file cannot be NULL" = !is.null(out_file))

    p = ggplot(df, aes(x=Count, y=fct_rev(chr), label=Count)) +
        geom_line(aes(size=Count*0.1), color="grey", alpha=0.75, show.legend=FALSE) +
        geom_point(aes(size=Count, color=ref_name), alpha=0.98, show.legend=TRUE) +
        labs(title=title_text,
             subtitle=subtitle_text,
             x="Counts", y="Chromosome",
             color="Reference", size="Counts") +
        geom_text_repel(size=3, min.segment.length=0.5) +
    
        # Map color to target_type with named vector; finer control of colors
        scale_color_manual(values=set_names(df$color, df$ref_name)) +  
        theme_bw() +
        theme(legend.position="right", plot.title.position="panel") +
        guides(size=FALSE)  # Remove size legend (Counts)

    ggsave(p, filename=out_file, device="png", units="in", width=w, height=h, dpi=300)
    cat("Wrote", out_file)
    return(list(path = out_file, plot = p))
}

#' @title Concatenate FASTA stats files
concat_fasta_stats <- function(file_paths, ref_names, ref_colors) {

    stopifnot(length(file_paths) == length(ref_names),
              length(ref_names) == length(ref_colors))
    
    df = map2(.x = file_paths, .y = ref_names, .progress = FALSE,
              .f = ~ fread(.x, data.table=FALSE) %>% mutate(ref_name = .y)) %>%
             bind_rows()
    
    df$color = ref_colors[match(df$ref_name, ref_names)]  # Add color column mapped by ref name
    
    # Keep FASTA stats for file with most bases in sum_len (i.e. genomic fasta) 
    df_stats = df %>% group_by(ref_name) %>%
        slice_max(order_by=sum_len, n=1) %>% ungroup() %>%
        mutate(chr="All\nContigs") %>% rename(Count = sum_n) %>%
        select(ref_name, color, Count, chr)
    return(df_stats)    
}

#' @title Get intersections
#'
#' @description Get binary intersection dataframe (1=Yes/Present, 0=No/Missing)
#'
#' @note df must contain columns: feat_id, ref_name
get_intersections <- function(df){
    df_binary = df %>% select(feat_id, ref_name) %>%
        mutate(value = 1) %>%
        pivot_wider(names_from = ref_name, values_from = value) %>%
        replace(is.na(.), 0) %>%
        column_to_rownames("feat_id")    
    return(df_binary)
}

#' @title Plot feature overlaps across multiple references
plot_common_features <- function(df, feat_type=NULL) {
    
    feat_type = arg_match(feat_type, c("gene_name", "gene_id", "transcript_id"), error_arg="feat_type")
    feature = str_split(feat_type, "_", simplify=TRUE)[1]
    feat_name = str_split(feat_type, "_", simplify=TRUE)[2]
    
    p = UpSetR::upset(
        df,
        keep.order = TRUE,  # Default: TRUE orders by set_cols, FALSE orders by set sizes
        empty.intersections = "on",  # Up to nintersects
        order.by = c("freq", "degree"),  # Freq orders high to low, degree orders by number of intersections
        show.numbers = "yes",  # Intersection sizes above bars
        number.angles = 0,  # For show.numbers
        set_size.show = TRUE,
        set_size.angles = 0,
    
        ## Layout parameters ##
        # c(I size title, I size tick labels, S size title, S size tick labels, S names, numbers above bars)
        text.scale = c(2, 1.7, 1.1, 1, 1.3, 1.1),  # Universal scale, or 6 individual sizes (see above)
        mainbar.y.label = paste0(str_to_title(feature),"s per Intersection\n(by ",feature," ",feat_name,")"),
        sets.x.label = paste0(str_to_title(feature),"s per Reference"),
        main.bar.color = "violetred4",
        sets.bar.color = "turquoise4",
        matrix.color = "slateblue4",
        shade.color = "grey60",
        mb.ratio = c(0.65, 0.35),  # default: c(0.7, 0.3); Ratio between matrix plot and main bar plot
        matrix.dot.alpha = 0.5,
        shade.alpha = 0.25,
        point.size = 3,  # default: 2.2; For matrix plot
        line.size = 0.7  # default: 0.7; For matrix plot
    )
    
    plot_name = paste0("Common ",feature,"s across ",ncol(df)," reference(s)")
    fig_name = paste0(make_clean_names(plot_name, replace=c(`\\(`="")), ".png")
    out_file = file.path(plot_path, fig_name)
    
    png(file=out_file, width=7, height=5.5, units="in", res=1000)
    print(p)
    dev.off()
    cat("Wrote", out_file)
    return(list(path = out_file, plot = p))
}

#' @title Summarize common features
#'
#' @description Get common features across multiple references
#'
#' @note df must contain columns: feat_id, ref_name, chr, strand, length
summarize_common_features <- function(df){
    
    total_refs = n_distinct(df$ref_name)
    
    df_common = df %>% select(feat_id, ref_name, chr, strand, length) %>%
        group_by(feat_id) %>%
        summarize(
            ref_names = paste(ref_name, collapse = "; "),
            across(all_of(c("chr","strand","length")),
                ~ ifelse(n_distinct(.x)==1,
                    paste(unique(.x), collapse = "; "),
                    paste(.x, collapse = "; ")), .names="{.col}s"),
              ref_counts = n(),
              across(all_of(c("chr","strand","length")), ~ n_distinct(.x), .names="{.col}_counts"),
              is_same_across_refs = case_when(
                  ref_counts == total_refs &
                  chr_counts == 1 &
                  strand_counts == 1 &
                  length_counts == 1
                  ~ "Yes", .default="No")
        ) %>% ungroup()
    return(df_common)
}

#' @title Create formatted Excel file
#'
#' @description Light-weight version of vp.create_excel_tab
#'
#' @note Add: Format multi-line headers (df_headers). Group columns (row_to_group). Cond. format (cnd_format_dict).
create_excel_tab <- function(df, tab_name="Sheet", out_file=NULL, prev_workbook=NULL, template_file=NULL,
                             plot_col=NULL, plot_row=NULL, freeze_cell=c(1,2), header_color="#F7FFBF") {
    
    # Use previous workbook, template, or create a new one from scratch
    if (!is.null(prev_workbook)) {
        wb = prev_workbook
    } else {
        if (!is.null(template_file)) {
            wb = loadWorkbook(template_file, isUnzipped=TRUE)
        } else {
            wb = createWorkbook() }}
    
    addWorksheet(wb, tab_name)
    writeData(wb, tab_name, df, withFilter=TRUE)
    
    # Set hyperlinks as formulas, and add as either column or row
    n_cols = ncol(df)
    if (!is.null(plot_col)) {
        n_cols = ncol(df)+1
        for (i in seq_along(plot_col)) {
            cell_value = plot_col[i]
            if (str_starts(cell_value, "=HYPERLINK")) {
                writeFormula(wb, tab_name, x=cell_value, startCol=n_cols, startRow=i)
            } else {
                writeData(wb, tab_name, x=cell_value, startCol=n_cols, startRow=i) }}}
    
    if (!is.null(plot_row)) {
        for (i in seq_along(plot_row)) {
            cell_value = plot_row[i]
            if (str_starts(cell_value, "=HYPERLINK")) {
                writeFormula(wb, tab_name, x=cell_value, startCol=i, startRow=nrow(df)+2)
            } else {
                writeData(wb, tab_name, x=cell_value, startCol=i, startRow=nrow(df)+2) }}}
    
    # Write Excel with header style and pane freeze (for scrolling)
    hs = createStyle(fgFill=header_color, halign="left", textDecoration="bold", border="TopBottomLeftRight")
    addStyle(wb, tab_name, style=hs, rows=1, cols=1:n_cols)
    freezePane(wb, sheet=tab_name, firstActiveCol=freeze_cell[1], firstActiveRow=freeze_cell[2])
    
    if (!is.null(out_file)) {
        saveWorkbook(wb, out_file, overwrite=TRUE)  # Write out final Excel file with all tabs
        cat("Wrote ***", out_file, " ***") }
    return(wb)
}

#' @title Plot genomic tracks
#'
#' @description Plot tracks for features across multiple references
#'
#' @note gggenomes uses hard-coded column names for dfs: seq_id, start, end, feat_id, length, strand
#'
#' @note Add: Plot multiple features together (including UTRs). Plot transcripts in region.
plot_track <- function(df, refs_to_plot=NULL, feats_to_plot=NULL, feat_type=NULL,
                       plot_path=NULL, color="basetheme::deepblue", show_track_info=FALSE,
                       w=5, h=5, v=1.25) {

    feat_types = c("gene_name", "gene_id", "transcript_id", "exon", "cds", "region")
    feat_type = arg_match(feat_type, feat_types, error_arg="feat_type")
    stopifnot("refs_to_plot cannot be NULL" = !is.null(refs_to_plot))
    stopifnot("feats_to_plot cannot be NULL" = !is.null(feats_to_plot))
    stopifnot("plot_path cannot be NULL" = !is.null(plot_path))
    
    if (feat_type != "region") {
        # Set filter column based on gene or transcript prefix in feats_to_plot
        if (any(str_detect(feats_to_plot, "^(ENST|NM_|NR_|XM_|XR_)"))){
            feat_col_name = sym("transcript_id")
        } else if (any(str_detect(feats_to_plot, "^(ENSG)"))){
            feat_col_name = sym("gene_id")
        } else {
            feat_col_name = sym("gene_name") }

        # Set feat_type differently for "exon" and "cds", since we want them on same seq_id (transcript) 
        if (feat_type %in% c("gene_name", "gene_id", "transcript_id")) {
            feat_type_sym = sym(feat_type)
        } else {
            feat_type_sym = feat_col_name }

        gene_track = df %>%
            filter({{feat_col_name}} %in% feats_to_plot, ref_name %in% refs_to_plot) %>%
            arrange({{feat_col_name}}, -length, {{feat_type_sym}}) %>%
            mutate(feat_id={{feat_type_sym}})
            
        if (feat_type %in% c("gene_name", "gene_id", "transcript_id")) {
            gene_track = gene_track %>%
                mutate(seq_id=paste0("[", ref_name, "] ", "[", feat_id, "] ", chr, ":", start, "-", end)) %>%
                mutate(seq_id=factor(seq_id, levels=.$seq_id))
        } else {
            gene_track = gene_track %>%
                mutate(seq_id=paste0("[", ref_name, "] ", "[", feat_id, "] ")) }
    
        feat_name = paste(feats_to_plot, collapse=", ")
        feature = str_split(feat_type, "_", simplify=TRUE)[1]
        plot_name = paste0(feat_name," ",feature,"s across ",length(refs_to_plot)," reference(s)")
        
    } else if (feat_type == "region") {
        split_region = str_split(feats_to_plot, "[:-]", simplify=TRUE)
        coi = split_region[1]
        soi = as.numeric(split_region[2])
        eoi = as.numeric(split_region[3])

        gene_track = df %>%
            filter(chr==coi & start>soi & end<eoi & ref_name %in% refs_to_plot) %>%
            mutate(feat_id=gene_name, seq_id=paste0("[", ref_name, "] "))
        
        plot_name = paste0("chr",feats_to_plot," across ",length(refs_to_plot)," reference(s)") }

    p = gggenomes(genes=gene_track) +
        geom_gene(aes(fill=feat_id), position="strand") +
        geom_seq_label(size=3, vjust=v) +
        scale_fill_paletteer_d(color)

    if (feat_type != "region") {
        p = p + labs(title=plot_name) + theme(legend.position="none")
        if (feat_type != "exon") {
            p = p + geom_gene_label(aes(label=paste0(length," bp")), angle=0, size=4, vjust=-1, hjust=0.5)
        } else {
            p = p + geom_gene_label(aes(label=exon_number), angle=0, size=4, vjust=-1, hjust=0.5) }
    } else {
        p = p + labs(title=plot_name, fill="Genes") +
            theme(legend.position="bottom") }

    if (show_track_info){
        IRdisplay::display(p %>% track_info)}  # df of counts for seqs, genes, and feats

    fig_name = paste0(make_clean_names(plot_name, replace=c(`\\(`="")), ".png")
    out_file = file.path(plot_path, fig_name)
    ggsave(p, filename = out_file, device="png", units="in", width=w, height=h, dpi=300)
    
    cat("Wrote", out_file)
    return(list(path = out_file, plot = p))
}
