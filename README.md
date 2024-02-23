# gocam_enrichment

This repository implements a step-centric enrichment analysis on Gene Ontology Causal Activity Models (GO-CAMs).

See **example.ipynb** under **notebooks**. We hope to integrate this tool soon into a Genome Alliance resource such that the backend database of GOCAM models can be maintained and updated. However, until then, users could follow example.ipynb to perform analyses with this tool.

Due to limitations on file size, **data** does not contain the reacto.owl file, but it can be downloaded at http://purl.obolibrary.org/obo/go/extensions/reacto.owl. This file is used to map reactome identifiers to UniProt accession numbers, which are used in GO-CAM ttl files. These can be accessed at https://github.com/geneontology/noctua-models or queried at  http://rdf.geneontology.org/blazegraph/sparql.

**dev** contains jupyter notebooks that were used to construct gene sets from GO-CAM pathways and python modules for performing enrichment.

notebooks pertaining to analysis are located in **notebooks**

