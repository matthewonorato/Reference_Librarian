import VinlandPy as vp

def get_urls(source, species, release, file_type, readme_file, release_type="latest"):
    """
    Get base URLs from Ensembl or RefSeq for the specified file type(s). Return filename:URL as dict.
        - **Ensembl:** Get URL for primary/top-level FASTA files (no haplotype/patches, unmasked, unique portion of Y)
        - **RefSeq:** Get URL for each primary chr (concat, unmasked, hard mask Y) or genomic.fna (soft-mask)

    *Ensembl*: Gene models are annotated on reference genome and include in silico predictions
        - Files are stored in separate folders based on file type (ex. fasta or GTF)
            - Human and mouse primary assembly files = other species top-level files
            - Choice of unmasked, soft-mask (sm; repetitive seqs in lower case), or hard-mask (rm; NNN)
            - Note: See README files in individual folders for more info (ex. README for GTF file contents)
        - More (only in some species): mysql (whole Dbs), variations/[vep|vcf|gvf], regulation, etc.
            - Note: VEP files can be large (ex. > 20Gb zipped) and should be installed as a local database/cache

    *RefSeq*: Gene models are annotated on mRNA sequence (stop codon may produce shorter protein) and include
              population-specific variants (may be harder to map genomic variants via SNP calling)
        - All files are stored in latest_assembly_versions/{release}
            - Note: URLs only for vertebrates (genomes/refseq/vertebrate_mammalian is hard-coded)
        - Note: If the release is unknown, set it to "" and click URL to find the latest assembly versions 
        - More (dbSNP, ClinVar): https://ftp.ncbi.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers

    :param source: Either "Ensembl" or "RefSeq"
    :type source: string
    :param species: Name of species formatted as "genus_species" (ex. homo_sapiens); Set to "" for Ensembl's vep
    :type species: string
    :param release: Ensembl or RefSeq release version (Ensembl ex. "111", RefSeq ex. "GCF_000001405.40_GRCh38.p14")
    :type release: string
    :param file_type: Name of file type (see config for Ensembl and RefSeq options) 
    :type file_type: string
    :param readme_file: README file name (ex. README_GCF_000001405.40-RS_2024_08 )
    :type readme_file: string
    :param release_type: *default: "latest"* | Either "latest" or "older" (needed for older releases)
    :type release_type: string
    :return: File names and URLs as key-value pairs
    :rtype: dict
    :raises ValueError: If source is not "Ensembl" or "RefSeq" (or any case-insensitive variant)
    """ 
    if source.lower()=="ensembl":
        print("~ Ensembl ~")
        print("All Downloads: https://useast.ensembl.org/info/data/ftp/index.html")
        species_stats = f"https://useast.ensembl.org/{species.capitalize()}/Info/Annotation"
        print(f"Species Stats: {species_stats}")

        if file_type=="all":
            file_type = "fasta_dna,fasta_cdna,fasta_cds,fasta_ncrna,fasta_pep,gtf,tsv,fasta_stats"
            if species.lower() not in ["homo_sapiens", "mus_musculus"]:
                file_type = file_type.replace("fasta_dna","fasta_dna_index")
        
    elif source.lower()=="refseq":
        species = species.capitalize()
        full_release = readme_file.replace("README_", "")
        refseq_base_url = "https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian"
        
        print("~ RefSeq ~")
        print("All Downloads: https://ftp.ncbi.nih.gov/genomes/refseq")
        species_stats = f"https://www.ncbi.nlm.nih.gov/genome/annotation_euk/{species}/{full_release}"
        print(f"Species Stats: {species_stats}")

        if file_type=="all":
            file_type = (
                "fasta_dna,fasta_mrna,fasta_prna,fasta_cds,fasta_pep,fasta_pep_cds,"
                "gtf,fasta_stats,readme,assembly_info"
            )
    else:
        raise ValueError(f"Unknown source: {source}")
    
    file_types = file_type.replace(" ","").split(",")

    url_dict = {}
    for file_type in file_types:
        if source.lower()=="ensembl":
            file_subtype = file_type.split("_")[-1]
            base_url = f"https://ftp.ensembl.org/pub/release-{release}/{file_type}/{species.lower()}/{file_subtype}"
        elif source.lower()=="refseq":
            if release_type=="latest":
                base_url = f"{refseq_base_url}/{species}/latest_assembly_versions/{release}"
            elif release_type=="older":
                base_url = f"{refseq_base_url}/{species}/annotation_releases/{full_release}"

            if file_type=="readme":
                base_url = f"{base_url}/{readme_file}"

        if file_type!="fasta_stats":
            url_dict[f"{file_type}"] = base_url
        else:
            url_dict[f"{file_type}"] = species_stats
    return url_dict

def get_primary_info(assembly_file, assembly_names=["Primary Assembly", "non-nuclear", "C57BL/6J"]):
    """
    Fill_in_ChatGPT

    *RefSeq only*: Get info for RefSeq primary assemblies, including UNL/UNP scaffolds. Return DataFrame.
    """
    if source.lower()=="ensembl":
        print(f"Ensembl genome for {species} is already primary assembly/top-level")
        return pd.DataFrame.from_dict({"Ensembl": ["This is intentionally empty"]})
    elif source.lower()=="refseq":
        df_chr = pd.read_csv(assembly_file, comment="#", header=None, sep="\t")

        # Keep chr_name[0],assembly_type[1],chr_number[2],refseq_accn[6],assembly_unit[7],chr_length[8],ucsc_name[9]
        col_dict = {0:"chr_name", 1: "assembly_type", 2:"chr_number",
                    6:"refseq_accn", 7:"assembly_unit", 8:"chr_length", 9:"ucsc_name"}
        df_chr = df_chr[df_chr[7].isin(assembly_names)][[0, 1, 2, 6, 7, 8, 9]].rename(columns=col_dict)

        df_chr["chr_file_name"] = "chr" + df_chr["chr_number"]
        df_chr["chr_file_name"] = df_chr["chr_file_name"].str.replace("chrna", "unplaced")
        df_chr.loc[df_chr["assembly_type"]=="unlocalized-scaffold", "chr_file_name"] = \
        df_chr["chr_file_name"] + ".unlocalized"

        out_file = assembly_file.replace("report.txt", "primary_info.csv")
        df_chr.to_csv(out_file, index=False)
        return df_chr

def download_files(source, url_dict):
    """
    Download files from Ensembl or Refseq based on config parameters. Return filename:filepath as dict.
        - Note: Additional files can be added later with `add_ref_to_master_file()`
    
    *Ensembl* primary/top-level.fa.gz are different sizes in dna and dna_index, but uncompressed they are the same
    
    *Ensembl* cDNA.fa is in the same sense orientation (5'→3') as *RefSeq* rna.fna (mRNA as it would be translated)

    Note: T2T's chm13v2.0_maskedY has chr1-22 renamed to 1-22 (https://github.com/marbl/CHM13)
    
    :param source: Either "Ensembl" or "RefSeq" (any case is fine)
    :type source: string
    :param url_dict: File names and URLs as key-value pairs
    :type url_dict: dict
    :return: File names and file paths as key-value pairs
    :rtype: dict
    """
    file_prefix = f"{species}.{genome_version}"
    file_full_prefix = f"{file_prefix}.{release}"
    
    file_dict = {}
    for file_type, base_url in url_dict.items():
        if source.lower()=="ensembl":
            # Genomic fasta of headers:sequences - primary assemblies only (human and mouse only)
            if file_type.lower()=="fasta_dna":
                gx_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.dna.primary_assembly.fa.gz"),
                                           out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = gx_file
                
            # Genomic fasta of headers:sequences - toplevel = primary (non-human and non-mouse species)
            if file_type.lower()=="fasta_dna_index": 
                gx_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.dna.toplevel.fa.gz"),
                                           out_path=PATH_TO_NEW_REFERENCE)
                file_dict["fasta_dna"] = gx_file  # Set to "fasta_dna" for easier downstream processing

            if file_type.lower()=="gtf":
                gtf_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.{release}.gtf.gz"),
                                            out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = gtf_file

            if file_type.lower()=="fasta_cdna":
                cdna_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.cdna.all.fa.gz"),
                                             out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = cdna_file

            if file_type.lower()=="fasta_cds": # Genomic CDS fasta - No UTR or intronic sequences
                cds_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.cds.all.fa.gz"),
                                            out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = cds_file

            if file_type.lower()=="fasta_ncrna": # *Transcript* sequences of ncRNA
                ncrna_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.ncrna.fa.gz"),
                                              out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = ncrna_file

            if file_type.lower()=="fasta_pep": # *Protein* translations of Ensembl genes
                aa_file = vp.download_file(os.path.join(base_url, f"{file_prefix}.pep.all.fa.gz"),
                                           out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = aa_file
            
            # RefSeq, Entrez, UniProt, ENA, and karyotype mapping files (ex. Ensembl and RefSeq transcript IDs)
            if file_type.lower()=="tsv":
                map_path = vp.create_dir(os.path.join(PATH_TO_NEW_REFERENCE, "mappings"), show_print=False) 
                file_dict[f"{file_type}"] = map_path

                refseq_file    = vp.download_file(os.path.join(base_url, f"{file_full_prefix}.refseq.tsv.gz"),
                                                  out_path=map_path)
                entrez_file    = vp.download_file(os.path.join(base_url, f"{file_full_prefix}.entrez.tsv.gz"),
                                                  out_path=map_path)
                uniprot_file   = vp.download_file(os.path.join(base_url, f"{file_full_prefix}.uniprot.tsv.gz"),
                                                  out_path=map_path)
                ena_file       = vp.download_file(os.path.join(base_url, f"{file_full_prefix}.ena.tsv.gz"),
                                                  out_path=map_path)
                karyotype_file = vp.download_file(os.path.join(base_url, f"{file_full_prefix}.karyotype.tsv.gz"),
                                                  out_path=map_path)
        elif source.lower()=="refseq":        
            if file_type.lower()=="fasta_dna":  # Genomic fasta with patches (sequence updates) and alt. chr.
                if "T2T" not in release:
                    gx_file = vp.download_file(os.path.join(base_url, f"{release}_genomic.fna.gz"),
                                               out_path=PATH_TO_NEW_REFERENCE)
                    file_dict[f"{file_type}"] = gx_file
            
            if file_type.lower()=="gtf":
                gtf_file = vp.download_file(os.path.join(base_url, f"{release}_genomic.gtf.gz"),
                                            out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = gtf_file
    
            if file_type.lower()=="fasta_mrna":  # Mature RNA fasta - exons + UTRs for mRNA, rRNAs, tRNAs
                mrna_file = vp.download_file(os.path.join(base_url, f"{release}_rna.fna.gz"),
                                             out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = mrna_file

            if file_type.lower()=="fasta_prna":  # Primary RNA fasta - exons + UTRs + introns for mRNA, rRNAs, tRNAs
                prna_file = vp.download_file(os.path.join(base_url, f"{release}_rna_from_genomic.fna.gz"),
                                             out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = prna_file
            
            if file_type.lower()=="fasta_cds":  # Genomic CDS fasta - No UTR or intronic sequences
                cds_file = vp.download_file(os.path.join(base_url, f"{release}_cds_from_genomic.fna.gz"),
                                            out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = cds_file
    
            if file_type.lower()=="fasta_pep":  # Protein fasta - Experimentally validated proteins
                aa_file = vp.download_file(os.path.join(base_url, f"{release}_protein.faa.gz"),
                                           out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = aa_file
    
            if file_type.lower()=="fasta_pep_cds":  # Proteomic CDS fasta - in silico translation from genomic DNA
                aa_cds_file = vp.download_file(os.path.join(base_url, f"{release}_translated_cds.faa.gz"),
                                               out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = aa_cds_file
            
            # Assembly metadata (ex. release date, accession:name pairs, features - Gene; mRNA; CDS; ncRNA)
            if file_type.lower()=="readme":
                readme_file = vp.download_file(base_url, out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = readme_file

            # Sequencing metadata and per chr. info (ex. chr. number and GenBank/RefSeq accessions)
            if file_type.lower()=="assembly_info":
                assembly_file = vp.download_file(os.path.join(base_url, f"{release}_assembly_report.txt"),
                                                 out_path=PATH_TO_NEW_REFERENCE)
                file_dict[f"{file_type}"] = assembly_file
            
    if source.lower()=="refseq" and "T2T" in release:       
        gx_file = os.path.join(PATH_TO_NEW_REFERENCE, f"{release}_genomic.fna")
        if os.path.isfile(gx_file):
            print(f"{gx_file} already exists")
            file_dict[f"fasta_dna"] = gx_file
        else:
            t2t_path = ("https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/"
                        "CHM13/assemblies/analysis_set/chm13v2.0_maskedY.fa.gz")
            gx_tmp_file = vp.download_file(t2t_path,
                                           out_path=PATH_TO_NEW_REFERENCE)

            df_chr = get_primary_info(assembly_file)
            sed_str = " ".join([f"-e 's/^>{chr_file_name}/>{chr_name}/'"
                                for chr_file_name, chr_name in zip(df_chr["chr_file_name"], df_chr["chr_name"])])
            
            sed_cmd = f"sed {sed_str} {gx_tmp_file} > {gx_file}"  # Replace chr_file_name with chr_name in new file
            vp.run_cmd(sed_cmd)
            os.remove(f"{gx_tmp_file}")
            file_dict[f"fasta_dna"] = gx_file      
    return file_dict

def update_chr_name(df_chr, chr_path):
    """
    Fill_in_ChatGPT
    
    *RefSeq only*: Convert RefSeq accession to chr number in fasta header using sed (ex. 'NC_052255.1' to '1')
    """
    chr_files = glob(os.path.join(chr_path, "*.fna"))
    chr_file_names = [f.split("/")[-1].rsplit(".", maxsplit=2)[0] for f in chr_files]  # Ex: (chr1.unlocalized).scaf.fna

    df_chr_file = pd.DataFrame.from_dict({"chr_file": chr_files, "chr_file_name": chr_file_names})
    df_chr_file = df_chr_file.merge(df_chr, how="left", on="chr_file_name")
    df_chr_file = df_chr_file[df_chr_file["refseq_accn"]!="na"]

    for chr_name in df_chr_file["chr_file_name"].unique():
        df_per_chr = df_chr_file[df_chr_file["chr_file_name"]==chr_name]
        chr_file = df_per_chr["chr_file"].iloc[0]  # File for new chr. name
        chr_ids = df_per_chr["refseq_accn"].str.replace(".", r"\.").tolist()  # Original name with escape char for sed
        chr_names = df_per_chr["chr_name"].tolist()  # New chr. name

        sed_str = " ".join([f"-e 's/^>{chr_id}/>{chr_name}/'" for chr_id, chr_name in zip(chr_ids, chr_names)])
        chr_new_file = chr_file.replace(".fna", ".renamed.fna")
        sed_cmd = f"sed {sed_str} {chr_file} > {chr_new_file}"  # sed substitute chr_id for chr_name in new file
        vp.run_cmd(sed_cmd)
    print("Converted headers from RefSeq Accn to chr name for RefSeq primary assemblies")

def mask_par_xy(base_url, df_chr, chr_path):
    """
    Fill_in_ChatGPT
    
    *RefSeq only*: Mask PAR region on Y chr, since it is a duplicate region of X chr. No return value.
       - Note: Ensembl uses unique region of Y chr (similar to masking)
    
    PAR regions: pseudoautosomal regions (PAR1 and PAR2) are short regions of homology b/w mammalian X and Y chrs
        - PAR behave like an autosome and recombine during meiosis
        - Thus genes in this region are inherited in an autosomal rather than a strictly sex-linked fashion

    Reasons to hard-mask
        - DNA-Seq: Reduces ambiguous alignments (+ mapping by up to 20%); Improves germline variant calling accuracy
        - RNA-Seq: Prevents misalignment of identical regions on chrY, which forces reads to align to chrX
    """
    if species.lower() in ["homo_sapiens", "mus_musculus"]:
        par_file = os.path.join(base_url,
                               f"{release}_assembly_structure/Primary_Assembly/pseudoautosomal_region/par_align.gff")
        f = requests.get(par_file)
        par_locs = [f.replace("Target=", "") for f in f.text.strip().split(";") if f.startswith("Target=")]

        par_bed = os.path.join(chr_path, "par.bed")
        chr_names = []
        with open(par_bed, "w") as out_file:
            for par_loc in par_locs:
                par_loc = par_loc.split(" ")
                chr_name = df_chr[df_chr["refseq_accn"]==par_loc[0]]["chr_name"].iloc[0]  # chr_name from X or Y
                par_loc[0] = chr_name
                par_loc[1] = str(int(par_loc[1])-1)
                out_file.write("\t".join(par_loc)+"\n")
                chr_names.append(chr_name)

        xy_file = os.path.join(chr_path, f"chr{chr_names[0]}.renamed.fna")  # Assumes PAR region on one chr. (all on Y)
        xy_mask_file = xy_file.replace(".fna", ".masked.fna")
        bedtools_path = "/mnt/disks/resources/software/bedtools2/bin/bedtools"
        mask_cmd = f"{bedtools_path} maskfasta -fi {xy_file} -bed {par_bed} -fo {xy_mask_file}"  # fi != fo
        vp.run_cmd(mask_cmd)

        os.system(f"cp {xy_mask_file} {xy_file}")
        os.remove(f"{xy_mask_file}")
        print(f"Used PAR gff/bed file to mask PAR region on {xy_file}")
        return

def download_primary_assemblies(url_dict, file_dict, delete_temp_folder=True, overwrite_file=True):
    """
    Fill_in_ChatGPT
    
    *RefSeq only*: Dowload primaries from NCBI, rename fasta headers, and mask PAR region of Y chr. This is done
    to create a FASTA file for primary chr. only (main chrs, MT, and UNL/UNP chrs). No return value.
        - Note: This will overwrite original genomic.fna with primary assemblies only
        - Note: Some species do not have an MT chr (like rat), so "MT" will not be present
    """
    if source.lower()=="ensembl":
        print(f"Ensembl genome for {species} is already primary assembly/top-level")
        return
    
    if source.lower()=="refseq" and "T2T" in release:
        print(f"RefSeq genome for {species} is already primary assembly/top-level")
        return

    gx_file = file_dict["fasta_dna"]
    if os.path.isfile(gx_file) and overwrite_file is False:  # File exists and we don't want to overwrite it
        print(f"{gx_file} already converted to primary assembly")
        return
    else:
        base_url = list(url_dict.values())[0]  # RefSeq URLs are all the same, so take 1st one
        out_path = vp.create_dir(os.path.join(PATH_TO_NEW_REFERENCE, "primary_assembly"), show_print=False)

        # Download files in specified directory (p_url, mt_url) without "robots" stopping download
        chr_types = ["assembled_chromosomes", "unlocalized_scaffolds", "unplaced_scaffolds"]  # URLs for each chr.
        for chr_type in chr_types:
            p_url = os.path.join(base_url, f"{release}_assembly_structure/Primary_Assembly/{chr_type}/FASTA/")
            dl_cmd = f'wget -r -l1 -e robots=off --no-parent --no-directories --directory-prefix="{out_path}" {p_url}'
            vp.run_cmd(dl_cmd)

        mt_url = os.path.join(base_url, f"{release}_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/")  
        dl_mt_cmd = f'wget -r -l1 -e robots=off --no-parent --no-directories --directory-prefix="{out_path}" {mt_url}'
        vp.run_cmd(dl_mt_cmd)

        gz_all_cmd = f"gunzip -r {out_path}"
        vp.run_cmd(gz_all_cmd)
        print(f"Unzipped and wrote primary assemblies for each chr. to {out_path}\n")

        assembly_file = file_dict["assembly_info"]
        df_chr = get_primary_info(assembly_file)
        update_chr_name(df_chr, out_path)  # Convert RefSeq-Accn to chr. name in fasta header
        mask_par_xy(base_url, df_chr, out_path)  # Human/Mouse: Mask PAR region on Y/X chr

        # Get all primary assembly files then concat in correct chr order (ex. chr1,chr2,...,chrX,chrY)
        chr_files = glob(os.path.join(out_path, "*.renamed.fna"))
        # Order: main chr (ex. 1-22), MT/X/Y chr, UNL on main chr, UNL on MT/X/Y chr, UNP chr
        sort_order = [os.path.join(out_path,f"chr{str(i)}.renamed.fna") for i in range(1,100)] +\
                     [os.path.join(out_path,f"chr{i}.renamed.fna") for i in ["MT", "X", "Y"]] +\
                     [os.path.join(out_path,f"chr{i}.unlocalized.scaf.renamed.fna") for i in range(1,100)] +\
                     [os.path.join(out_path,f"chr{i}.unlocalized.scaf.renamed.fna") for i in ["MT", "X", "Y"]] +\
                     [os.path.join(out_path,"unplaced.scaf.renamed.fna")]
        chr_sort_files = [f for f in sort_order if f in chr_files]
        assert len(chr_files) == len(chr_sort_files)

        cat_cmd = f"cat {' '.join(chr_sort_files)} > {gx_file}"
        vp.run_cmd(cat_cmd)
        print(f"Wrote {gx_file}")

        if delete_temp_folder:
            rm_cmd = f"rm -r {out_path}"
            vp.run_cmd(rm_cmd)
            print(f"Deleted {out_path}")
        return

def load_gtf(gtf_file):
    """
    Load in GTF file as polars df then convert to pandas df. Return df.

    Genomic GTF files contain ALL gene, transcript, and exon annotations (but no sequences)
        - More: https://www.gencodegenes.org/pages/faq.html or https://genome.ucsc.edu/FAQ/FAQgenes.html
        - Note: For Ensembl and RefSeq, missing gene names are due to "novel genes" without names

    Additional GTF/GFF tools:
        - (good features - get gene lengths, find retained introns, etc) https://github.com/RacconC/gtftools
        - (comprehensive) https://github.com/NBISweden/AGAT

    :param file_dict: File names and file paths as key-value pairs
    :type file_dict: dict
    :return: DataFrame containing all GTF info (features, sources, coords, strand, names, etc.)
    :rtype: Pandas DataFrame
    :raises TypeError: If the file type is not supported
    """
    if not gtf_file.endswith(".gtf"):
        raise TypeError(f"{gtf_file} does not have a .gtf extension. Please provide another file.")

    # Convert polars df to pandas df; Rename seqname to chr (since we'll use chr in all subsequent files)
    df = gtfparse.read_gtf(gtf_file).to_pandas(use_pyarrow_extension_array=True).rename(columns={"seqname": "chr"})

    if "gene_name" not in df.columns:  # For RefSeq
        df = df.rename(columns={"gene": "gene_name"})
    return df

def get_anno_from_gtf(df_gtf, file_dict, ann_path, feature, cols, overwrite_file=False):
    """
    Fill_in_ChatGPT
    """
    feat_name = feature.replace("five_prime","5p").replace("three_prime","3p").lower()

    if source.lower()=="ensembl":
        bed_file = os.path.join(ann_path, f"{species}.{genome_version}.{release}.{feat_name}.bed")
    elif source.lower()=="refseq":
        bed_file = os.path.join(ann_path, f"{release}.{feat_name}.bed")

    out_file = bed_file.replace(".bed", ".sorted.tsv")
    if os.path.isfile(out_file) and overwrite_file is False:  # File exists and we don't want to overwrite it
        print(f"{out_file} already exists")
        return out_file
        
    df = df_gtf[df_gtf["feature"]==feature][cols] 
    df["score"] = 0  # Must be an integer b/w 0 and 1000
    
    if source.lower()=="refseq":  # Get chr name from RefSeq Accn (ex. NC_000001.11 = chr1)
        assembly_file = file_dict["assembly_info"]
        df_chr = get_primary_info(assembly_file)
        chr_dict = dict(zip(df_chr["refseq_accn"], df_chr["chr_name"]))
        df["chr"] = df["chr"].map(chr_dict)
        
    if feature!="gene":  # Add gene types to Tx and exons
        gene_types = df_gtf[df_gtf["feature"]=="gene"][["gene_id", "gene_biotype"]]
        gene_types_dict = dict(zip(gene_types["gene_id"], gene_types["gene_biotype"]))
        df["gene_biotype"] = df["gene_id"].map(gene_types_dict)

    df = df.dropna(subset=["chr"])  # Remove genes on non-primaries (genes on fix-patches or alt-scaffolds)
    df.to_csv(bed_file, sep="\t", header=False, index=False)

    # Remove genes on non-primaries (genes on fix-patches or alt-scaffolds)
    filter_cmd = f"awk '$1 ~ /^(chr)?([0-9]+|X|Y)$/' {bed_file}"
    
    # Sort BED file by chr (1-22, X, Y, MT) then gene start (numeric order) then gene end (numeric order)
    # -k1,1fV: sort chr. (key 1)  by ignorning case (f: MT after chr) and using natural sort (V: chr2 before chr10)
    sort_cmd = f"sort -k1,1fV -k2,2n -k3,3n > {out_file}"
    vp.run_cmd(f"{filter_cmd} | {sort_cmd}")

    header_col = "\\t".join([c for c in cols])
    sed_cmd = f"sed -i '1i\\{header_col}' {out_file}"
    vp.run_cmd(sed_cmd)

    print(f"Wrote {out_file}")
    os.remove(f"{bed_file}")
    return out_file

def create_anno_files(file_dict, ann_type, overwrite_files=False):
    """
    Fill_in_ChatGPT

    Get feature annotations based on config parameters. Annotations are extracted from GTF file and then sorted 
    by chromosome, similar to a BED file. Return filename:filepath as dict.

    - Note: Annotations for primary assemblies only (NC in Refseq, 1/X/MT in Ensembl)
    - *Ensembl only*: Only Ensembl has GTF annotations for 5'-UTRs and 3-UTRs
    - Anno vs. BED Files: one row per exon/contiguous region (long) vs. one row per feature (wide)
    """
    out_path = os.path.join(PATH_TO_NEW_REFERENCE, "annotations")
    vp.create_dir(out_path, show_print=False)
     
    df_gtf = load_gtf(file_dict["gtf"])

    if ann_type=="all":
        if source.lower()=="ensembl":
            ann_type = "gene,transcript,exon,cds,5p_utr,3p_utr"
        elif source.lower()=="refseq":
            ann_type = "gene,transcript,exon,cds"
    
    ann_types = ann_type.replace(" ","").split(",")

    ann_dict = {}
    for ann_type in ann_types:
        if ann_type.lower()=="gene":
            feat_to_use = ann_type
            cols = ["chr","start","end","gene_id","score","strand",
                    "gene_name","gene_biotype","source","frame"]
            
        if ann_type.lower()=="transcript":
            feat_to_use = ann_type
            cols = ["chr","start","end","transcript_id","score","strand",
                    "gene_id","gene_name","gene_biotype","transcript_biotype","source","frame"]
            
        if ann_type.lower()=="exon":
            feat_to_use = ann_type
            cols = ["chr","start","end","transcript_id","score","strand","exon_number",
                    "gene_id","gene_name","gene_biotype","transcript_biotype","source","frame"]
            
        if ann_type.lower()=="cds":
            feat_to_use = ann_type.upper()
            cols = ["chr","start","end","transcript_id","score","strand",
                    "gene_id","gene_name","gene_biotype","transcript_biotype","source","frame"]
            
        if "utr" in ann_type.lower() and source.lower() == "ensembl":
            cols = ["chr","start","end","transcript_id","score","strand",
                    "gene_id","gene_name","gene_biotype","transcript_biotype","source","frame"]
            
            if ann_type.lower()=="5p_utr":
                feat_to_use = "five_prime_utr"
            if ann_type.lower()=="3p_utr":
                feat_to_use = "three_prime_utr"

        ann_file = get_anno_from_gtf(
            df_gtf, file_dict, out_path,
            feature=feat_to_use, cols=cols,
            overwrite_file=overwrite_files
        )
        ann_dict[f"ann_{ann_type}"] = ann_file
    return ann_dict

def create_bed_files(file_dict, ann_dict, bed_type, overwrite_files=False):
    """
    Fill_in_ChatGPT

    Create BED files based on config parameters. In each BED file, feature coordinates are extracted from
    GTF file and then sorted by chromosome. Return filename:filepath as dict.

    - Note: BED files include primary assemblies only (NC in Refseq, 1/X/MT in Ensembl)
    - Rule: In BED file, start is always less than end regardless of strand (+ or -)
    
    - Note: IDs are preferable to names (ex. Ensembl has 60,000 unique IDs to 40,000 unique names)
        - For gene BED files, Ensembl uses gene ID (ENSG*) while RefSeq only has gene name (no IDs)

    bedparse tutorial & additional commands:
        - https://bedparse.readthedocs.io/en/latest/Tutorial.html
            - promoters: get region +/-NT from 0 for each transcript
            - validateFormat: check BED format or fix delims
            - convertChr: "1" to "chr1"
            - filter/join: by name (column 4)
    
    Good review of BED file:
        - https://samtools.github.io/hts-specs/BEDv1.pdf
        - BED12 format: https://bedtools.readthedocs.io/en/latest/content/general-usage.html

    - QC: Gene BED files do not have any missing names or duplicate names (true for Ensembl and RefSeq)
    - QC: bedparse bed12tobed6 creates same file as exons.sorted.tsv (checked + and - strands)
    """
    bed_path = os.path.join(PATH_TO_NEW_REFERENCE, "bed_files")
    vp.create_dir(bed_path, show_print=False)

    if source.lower()=="ensembl":
        file_full_prefix = f"{species}.{genome_version}.{release}"
    elif source.lower()=="refseq":
        file_full_prefix = f"{release}"

    # BED12 used as input to create other BED files
    tx_bed_file = os.path.join(bed_path, f"{file_full_prefix}.transcript.sorted.bed12")

    if bed_type=="all":
        if source.lower()=="ensembl":
            bed_type = "gene,transcript,exon,cds,5p_utr,3p_utr"
        elif source.lower()=="refseq":
            bed_type = "gene,transcript,exon,cds"

    bed_dict = {}
    bed_types = bed_type.replace(" ","").split(",")
    bedparse_path = "/mnt/disks/resources/software/miniconda3/envs/Py-3.12/bin/bedparse"
    
    for bed_type in bed_types:
        if bed_type.lower()=="transcript":
            out_file = tx_bed_file
        else:
            out_file = os.path.join(bed_path, f"{file_full_prefix}.{bed_type}.sorted.bed")
    
        if os.path.isfile(out_file) and overwrite_files is False:  # File exists and we don't want to overwrite it
            print(f"{out_file} already exists")
            bed_dict[f"bed_{bed_type}"] = out_file
        else:
            tmp_file = os.path.join(bed_path, f"{file_full_prefix}.{bed_type}.bed")

            if bed_type.lower()=="gene":
                gene_ann_file = ann_dict["ann_gene"]
                cols=["chr","start","end","gene_id","score","strand"]
                df_gene_ann = pd.read_csv(gene_ann_file, sep="\t")[cols]
                df_gene_ann["start"] = df_gene_ann["start"]-1  # FASTA to BED: 1-indexed to 0-indexed
                df_gene_ann.to_csv(out_file, sep="\t", header=False, index=False)
            if bed_type.lower()=="transcript":
                bed_cmd = f"{bedparse_path} gtf2bed {gtf_file} > {tmp_file}"
            if bed_type.lower()=="exon":
                bed_cmd = f"{bedparse_path} bed12tobed6 --whichExon all --appendExN {tx_bed_file} > {tmp_file}"
            if bed_type.lower()=="cds":
                bed_cmd = f"{bedparse_path} cds {tx_bed_file} > {tmp_file}"
            if bed_type.lower()=="5p_utr":
                bed_cmd = f"{bedparse_path} 5pUTR {tx_bed_file} > {tmp_file}"
            if bed_type.lower()=="3p_utr":
                bed_cmd = f"{bedparse_path} 3pUTR {tx_bed_file} > {tmp_file}"    

            if bed_type.lower()!="gene":
                vp.run_cmd(bed_cmd)

            # Since RefSeq stores Accn instead of chr in FASTA and GTF files, replace it (ex. NC_000001.11 = 1)
            if source.lower()=="refseq" and bed_type.lower()=="transcript":
                assembly_file = file_dict["assembly_info"]
                df_chr = get_primary_info(assembly_file)
                sed_str = " ".join([f"-e 's/^{chr_id}/{chr_name}/'"
                                    for chr_id, chr_name in zip(df_chr["refseq_accn"], df_chr["chr_name"])])
                
                sed_cmd = f"sed {sed_str} {tmp_file}"  # sed substitute chr_id for chr_name in new file
                filter_cmd = f"awk '$1 ~ /^(chr)?([0-9]+|X|Y)$/'"
                sort_cmd = f"sort -k1,1fV -k2,2n -k3,3n > {out_file}"
                vp.run_cmd(f"{sed_cmd} | {filter_cmd} | {sort_cmd}") 
            elif bed_type.lower()!="gene":   
                filter_cmd = f"awk '$1 ~ /^(chr)?([0-9]+|X|Y)$/' {tmp_file}"
                sort_cmd = f"sort -k1,1fV -k2,2n -k3,3n > {out_file}"
                vp.run_cmd(f"{filter_cmd} | {sort_cmd}")

            print(f"Wrote {out_file}")     
            bed_dict[f"bed_{bed_type}"] = out_file
            
            if os.path.isfile(tmp_file):
                os.remove(f"{tmp_file}")
    return bed_dict

def get_fasta_from_bed(file_dict, bed_dict, feature, n_threads=1, overwrite_file=False):
    """
    Fill_in_ChatGPT
    
    Extract feature from BED file in FASTA format. Return filename:filepath as dict.
        - Feature options: gene, transcript, exon, cds, 5p_utr, 3p_utr
        - BED: 0-based for start, 1-based for end (FASTA: chr1 1000 to 2000, BED: chr1 999 2000)
    
    Alt. Command: `bedtools getfasta` does not update FASTA coordinates from BED file
        - {bedtools_path} getfasta -name -s -split -fi {fasta_file} -bed {bed_file} -fo {new_fasta_file}
            - ID/name as header (-name), feats on - strand rev. comp. (-s), & exons/blocks concat (-split)
            - Note: If input FASTA file is RNA, add "-rna" option to command

    Note: Write mRNA/pRNA FASTA files only when needed since they are >10GB (ex. InSilico_OligoScreen)
    """
    fasta_file = file_dict["fasta_dna"]
    bed_file = bed_dict[f"bed_{feature}"]

    fasta_ext = fasta_file.split(".")[-1]
    new_fasta_name = bed_file.split("/")[-1].replace("sorted.bed",f"{fasta_ext}")
    new_fasta_file = os.path.join(PATH_TO_NEW_REFERENCE, f"{new_fasta_name}")

    if os.path.isfile(new_fasta_file) and overwrite_file is False:
        print(f"{new_fasta_file} already exists")
        return {f"fasta_{feature}": new_fasta_file}
    else:
        fasta_cmd = (
            f"gzip --stdout {fasta_file} | "  # zip and pass to stnd. out to retain order in BED file
            f"/mnt/disks/resources/software/seqkit/seqkit subseq "
            f"--bed {bed_file} "
            f"--out-file {new_fasta_file} "
            f"--line-width 0 "  # Create single-line sequences
            f"--threads {n_threads} "
            f"--seq-type auto "  # sequence type (dna|rna|protein|unlimit|auto); auto by first sequence
            f"--update-faidx"
        )
        vp.run_cmd(f"{fasta_cmd}")
    
        print(f"Wrote {new_fasta_file}")
        return {f"fasta_{feature}": new_fasta_file}

def get_fasta_stats(file_dict, n_threads=1, overwrite_file=False):
    """
    Fill_in_ChatGPT

    https://bioinf.shenwei.me/seqkit/usage/#stats
    """
    if source.lower()=="ensembl":
        out_file = os.path.join(PATH_TO_NEW_REFERENCE, f"{species}.{genome_version}.fasta_stats.tsv")
    elif source.lower()=="refseq":
        out_file = os.path.join(PATH_TO_NEW_REFERENCE, f"{release}.fasta_stats.tsv")
    
    if os.path.isfile(out_file) and overwrite_file is False:
        print(f"{out_file} already exists")
        return out_file
    else:
        fasta_list = []
        for file_type, file_path in file_dict.items():
            if "fasta" in file_type:
                fasta_list.append(file_path)
        
        stats_cmd = (
            f"/mnt/disks/resources/software/seqkit/seqkit stats "
            f"--all "
            f"--tabular "
            f"--threads {n_threads} "
            f"{" ".join(fasta_list)} "
            f"> {out_file}"
        )
        vp.run_cmd(stats_cmd)
        print(f"Wrote FASTA stats to {out_file}")
        return out_file

def create_genomic_index(fasta_file, overwrite_file=False):
    """
    Index FASTA file for fast access to data (ex. for IGV). Return path to index.

    :param fasta_file: Path to FASTA file
    :type fasta_file: string
    :param overwrite_file: *default: False* | Whether to overwrite pre-existing index file or not
    :type overwrite_file: Boolean
    :return: Path to FASTA index
    :rtype: string    
    """
    out_file = fasta_file + ".fai"
    
    if os.path.isfile(out_file) and overwrite_file is False:  # File exists and we don't want to overwrite it
        print(f"{out_file} already exists")
        return out_file
    else:  # File doesn't exist, or we want to overwrite it
        idx_cmd = f"/mnt/disks/resources/software/samtools-1.20/samtools faidx {fasta_file} -o {out_file}"
        vp.run_cmd(idx_cmd)
        print(f"Wrote {out_file}")
        return out_file

def create_bowtie_index(fasta_file, out_path, n_threads=1, overwrite_file=False):
    """
    Fill_in_ChatGPT

    https://bowtie-bio.sourceforge.net/manual.shtml
        - bowtie-build builds a Bowtie index from a set of DNA sequences (6 files, [1-6].ebwt)
        - These files together constitute the index: they are all that is needed to align reads to that reference
        - The original sequence files are no longer used by Bowtie once the index is built
        - By default, bowtie-build auto-searches for settings that yield the best running time w/o exhausting memory
    """
    vp.create_dir(out_path, show_print=False)

    # Files already in directory and we don't want to overwrite them
    files_in_dir = [f for f in os.listdir(out_path) if os.path.isfile(os.path.join(out_path, f))]
    if files_in_dir and overwrite_file is False:
        print(f"{out_path} already contains files")
        return out_path
    else:
        idx_name = fasta_file.rsplit(".", maxsplit=1)[0].split("/")[-1]
        idx_cmd = (
            f"/mnt/disks/resources/software/bowtie-1.3.1-linux-x86_64/bowtie-build "
            f"--threads {n_threads} "
            f"-f {fasta_file} "
            f"{out_path}/{idx_name}"
        )
        vp.run_cmd(idx_cmd)
        print(f"Wrote bowtie index to {out_path}")
        return out_path

def create_star_index(fasta_file, gtf_file, out_path, n_threads=1, overwrite_file=False):
    """
    Fill_in_ChatGPT
    """
    vp.create_dir(out_path, show_print=False)
    
    files_in_dir = [f for f in os.listdir(out_path) if os.path.isfile(os.path.join(out_path, f))]
    if files_in_dir and overwrite_file is False:
        print(f"{out_path} already contains files")
        return out_path
    else:
        idx_cmd = (
            f"/mnt/disks/resources/software/STAR-2.7.11b/source/STAR "
            f"--runMode genomeGenerate "
            f"--runThreadN {n_threads} "
            f"--genomeDir {out_path} "
            f"--genomeFasta_files {fasta_file} "  #  Should be plain text FASTA files; cannot be zipped?
            f"--sjdbGTFfile {gtf_file} "
            f"--sjdbOverhang 100"  # Ideal value (reads of varying length): {int(max(read_lengths)-1)}
        )
        vp.run_cmd(idx_cmd)
        print(f"Wrote STAR index to {out_path}")
        return out_path

def create_rsem_index(fasta_file, gtf_file, out_path, n_threads=1, overwrite_file=False):
    """
    Fill_in_ChatGPT
    """
    vp.create_dir(out_path, show_print=False)
    
    files_in_dir = [f for f in os.listdir(out_path) if os.path.isfile(os.path.join(out_path, f))]
    if files_in_dir and overwrite_file is False:
        print(f"{out_path} already contains files")
        return out_path
    else:
        idx_name = fasta_file.rsplit(".", maxsplit=1)[0].split("/")[-1]
        idx_cmd = (
            f"/mnt/disks/resources/software/RSEM-1.3.3/rsem-prepare-reference "
            f"--gtf {gtf_file} "
            f"{fasta_file} "
            f"{out_path}/{idx_name}"
        )
        vp.run_cmd(idx_cmd)
        print(f"Wrote RSEM index to {out_path}")
        return out_path

def create_salmon_index(fasta_file, out_path, n_threads=1, overwrite_file=False):
    """
    Fill_in_ChatGPT

    Salmon can use either Ensembl's cDNA.fa or RefSeq's rna.fna as both are in the sense orientation
    """
    vp.create_dir(out_path, show_print=False)
    
    files_in_dir = [f for f in os.listdir(out_path) if os.path.isfile(os.path.join(out_path, f))]
    if files_in_dir and overwrite_file is False:
        print(f"{out_path} already contains files")
        return out_path
    else:
        idx_cmd = (
            f"/mnt/disks/resources/software/miniconda3/envs/Py-3.12/bin/salmon index "
            f"--transcripts {fasta_file} "
            f"--index {out_path} "
            f"--threads {n_threads}"
        )
        vp.run_cmd(idx_cmd)
        print(f"Wrote Salmon index to {out_path}")
        return out_path

def create_indexes(file_dict, idx_type, overwrite_files=False):
    """
    Fill_in_ChatGPT

    Create indexes based on config parameters. Indexes allow us to quickly retrieve sequence data for
    specific regions or IDs, rather than having to scan the entire file each time (and load into memory).
    Return index_name:index_path as dict.
        - Example: Aligners like Bowtie can efficiently map reads to reference genomes using indexed sequences
        - Note: Chromosome names must match (FASTA chr1 = GTF chr1). If needed, use bedparse's 'convertChr'.
    """
    gx_file = file_dict["fasta_dna"]
    gtf_file = file_dict["gtf"]
    
    if source.lower()=="ensembl":
        tx_file = file_dict["fasta_cdna"]
    elif source.lower()=="refseq":
        tx_file = file_dict["fasta_mrna"]

    if idx_type=="all":
        idx_type = "fasta,bowtie_dna,bowtie_rna,star,rsem,salmon"
    
    idx_types = idx_type.replace(" ","").split(",")

    idx_dict = {}
    for idx_type in idx_types:
        if idx_type.lower()=="fasta":
            for file_type, file_path in file_dict.items():
                if "fasta" in file_type:
                    idx_file = create_genomic_index(file_path, overwrite_file=overwrite_files)
                    idx_dict[file_type.replace("fasta","idx")] = idx_file
    
        if idx_type.lower()=="bowtie_dna":
            bwt_dna_path = create_bowtie_index(
                fasta_file = gx_file,
                out_path = os.path.join(PATH_TO_NEW_REFERENCE, "bowtie_index", "dna"),
                n_threads = n_threads,
                overwrite_file = overwrite_files
            )
            idx_dict["idx_bwt_dna"] = bwt_dna_path
        
        if idx_type.lower()=="bowtie_rna":
            bwt_rna_path = create_bowtie_index(
                fasta_file = tx_file,
                out_path = os.path.join(PATH_TO_NEW_REFERENCE, "bowtie_index", "rna"),
                n_threads = n_threads,
                overwrite_file = overwrite_files
            )
            idx_dict["idx_bwt_rna"] = bwt_rna_path
        
        if idx_type.lower()=="star":
            star_path = create_star_index(
                fasta_file = gx_file,
                gtf_file = gtf_file,
                out_path = os.path.join(PATH_TO_NEW_REFERENCE, "star_index"),
                n_threads = n_threads,
                overwrite_file = overwrite_files
            )
            idx_dict["idx_star"] = star_path
    
        if idx_type.lower()=="rsem":
            rsem_path = create_rsem_index(
                fasta_file = gx_file,
                gtf_file = gtf_file,
                out_path = os.path.join(PATH_TO_NEW_REFERENCE, "rsem_index"),
                n_threads = n_threads,
                overwrite_file = overwrite_files
            )
            idx_dict["idx_rsem"] = rsem_path
    
        if idx_type.lower()=="salmon":
            salmon_path = create_salmon_index(
                fasta_file = tx_file,
                out_path = os.path.join(PATH_TO_NEW_REFERENCE, "salmon_index"),
                n_threads = n_threads,
                overwrite_file = overwrite_files
            )
            idx_dict["idx_salmon"] = salmon_path
    return idx_dict

def check_sums(in_file):
    """
    Fill_in_ChatGPT
    
    Update to accomodate Ensembl and RefSeq file validation checks (sum, etc)
    
    Ensembl's checksums uses `sum` for simple integrity checks (not MD5 hashes) on *zipped* files (fa.gz)
        - Outputs checksum_value and n_blocks (file size in GB = n_blocks * 512 bytes/block * GB/1e9 bytes)
    """
    return vp.run_cmd(f"md5sum {in_file}").stdout

def add_ref_to_master_file(ref_dict, species, new_file_url=None, new_ref_key=None):
    """
    Fill_in_ChatGPT

    Add files for species based on config parameters to master reference file. Return JSON file.
    """
    if not os.path.exists(PATH_TO_REFERENCES):
        with open(PATH_TO_REFERENCES, "w") as infile:
            json.dump(dict(), infile, indent=4)
    
    with open(PATH_TO_REFERENCES) as r:
        ref = json.load(r)
    
    # Save original, just in case
    og_file = PATH_TO_REFERENCES.replace(".json", "_original.json")
    with open(og_file, "w") as infile:
        json.dump(ref, infile, indent=4)
    
    species = species.capitalize()
    file_full_prefix = f"{species}.{genome_version}.{release}"
    
    if source.lower()=="ensembl":
        build_version = f"ensembl.{release}"
    elif source.lower()=="refseq":
        build_version = f"refseq.{release.rsplit("_", maxsplit=1)[0]}"
    
    new_ref_dict = {k:v for k,v in ref_dict.items() if ".fai" not in v}  # Remove .fai files from final dict
    
    if species not in ref.keys():  # Add new species
        ref[f"{species}"] = {f"{genome_version}": {f"{build_version}": new_ref_dict}}
        print(f"Added new species: {species}\n")
        
    elif genome_version not in ref[f"{species}"].keys():  # Add new genome
        ref[f"{species}"][f"{genome_version}"] = {f"{build_version}": new_ref_dict}
        print(f"Added new genome: {species}.{genome_version}\n")
        
    elif build_version not in ref[f"{species}"][f"{genome_version}"].keys():  # Add new build version
        ref[f"{species}"][f"{genome_version}"][f"{build_version}"] = new_ref_dict
        print(f"Added new build version: {file_full_prefix}\n")
        
    else:
        print(f"{file_full_prefix} already in config. No changes were made.\n")
    
    if new_file_url:  # Add new file
        new_file = vp.download_file(new_file_url, out_path=PATH_TO_NEW_REFERENCE)
        new_ref_dict = {new_ref_key: new_file}
        for k, v in new_ref_dict.items():
            if k not in ref[f"{species}"][f"{genome_version}"][f"{build_version}"].keys():
                ref[f"{species}"][f"{genome_version}"][f"{build_version}"].update({k: v})
                print(f"Added new file for: {file_full_prefix}\n")
        
    with open(PATH_TO_REFERENCES, "w") as in_file:
        json.dump(ref, in_file, indent=4)
    
    with open(PATH_TO_REFERENCES) as r:
        ref = json.load(r)
    return ref
