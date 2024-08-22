# Step Enrichment

Markarian N, Van Auken KM, Ebert D, Sternberg PW (2024) Enrichment on steps, not genes, improves inference of differentially expressed pathways. PLOS Computational Biology 20(3): e1011968. https://doi.org/10.1371/journal.pcbi.1011968

Associated notebooks can be found at https://github.com/nmarkari/gocam_enrichment

**This package also requires installation of R. Download the R package BiasedUrn package: https://cran.r-project.org/web/packages/BiasedUrn/BiasedUrn.pdf. This code was tested with R version 4.2.2 and BiasedUrn_2.0.11.**

Post-publication, we extended our code to work with KEGG data to aid other developers who wanted to use our tool, where KEGG reactions are used analagously to "sets" from our paper. Bulk downloads from KEGG require a subscription, so we are not redistributing their data here. You may download it through KEGGs FTP site if you do not already have access to the data via https://www.kegg.jp/kegg/download/. To use our kegg enrichment function, please create a "data/kegg/" directory and include the following files:

- ko-to-reaction.map
- module_name.txt
- pathway_name.txt
- reaction-to-ko.map
- reaction-to-module.map
- reaction-to-pathway.map
- reaction_enzyme.txt

Usage:\
\
enrich(filename, input_type, method = 'set', return_all = False, FDR=.05,fpath= '', display_gene_symbol = True, display_input = False, 
                   kegg = False, enrich_against = 'module'):\
\
-     method: set (unweighted enrichment on sets as described in our paper), ncHGT (weighted enrichment as described in our paper), or standard
-     valid gene IDs for GOCAMs/Reactome: 'Gene Symbol', 'ENSEMBL ID', 'uniprot'
-     valid gene IDs for kegg: 'ko', 'enzyme' (EC number)
- -       if enriching with kegg, set kegg = True.
- -       enrich against can be 'module' or 'pathways' for kegg
-     return_all: 
- -         if false, only returns the dataframe displaying results. \
- -         if true: returns (gene_list, filtered_out_genes, filtered_list, setID2members_input_uni, setID2members_input, df_display)
- -         return_all = True is not just for debugging. User may want to know which of their input genes were filtered out as well as how
- -         the IDs were mapped, as uniprot IDs can sometimes map to more than one HGNC gene symbol
-     FDR: false discovery rate (Benjamini Hochberg correction)
-     display_gene_symbol: if true, display HGNC symbols on output regardless of input ID type
-     display_input: overrides display_gene_symbol if true. If False, either displays the backend id type (UniProtKB IDs for REACTOME/GOCAMs, reactions for KEGG) or the default of HGNC gene names"""\
