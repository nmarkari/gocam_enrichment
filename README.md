# gocam_enrichment

**data** contains files pertaining to the contstruction of the reference gene sets for GO-CAM pathways, while **test_data** contains datasets that the tool was tested on. 

Due to limitations on file size, **data** does not contain the reacto.owl file, but it can be downloaded at http://purl.obolibrary.org/obo/go/extensions/reacto.owl. This file is used to map reactome identifiers to UniProt accession numbers, which are used in GO-CAM ttl files. These can be accessed at https://github.com/geneontology/noctua-models or queried at  http://rdf.geneontology.org/blazegraph/sparql.

**dev** contains jupyter notebooks that were used to construct gene sets from GO-CAM pathways and python modules for performing enrichment.
