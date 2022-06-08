=======================================
Explore the pangenome - PanTools manual
=======================================

.. container:: wy-grid-for-nav

   .. container:: wy-side-scroll

      .. container:: wy-side-nav-search

         `PanTools manual <..>`__

         .. container::

      .. container:: wy-menu wy-menu-vertical

         -  `Home <..>`__

         -  `Install <../install/>`__

         -  `Construct pangenome <../construct/>`__

         -  `Pangenome characterization <../characterize/>`__

         -  `Phylogeny <../phylogeny/>`__

         -  `Sequence alignments <../msa/>`__

         -  `Explore the pangenome <./>`__

            -  `Gene locations <#gene-locations>`__

               -  `Required arguments <#required-arguments>`__
               -  `Optional arguments <#optional-arguments>`__
               -  `Example command <#example-command>`__
               -  `Output files <#output-files>`__

            -  `Find genes <#find-genes>`__

               -  `Find genes by name <#find-genes-by-name>`__

                  -  `Required arguments <#required-arguments_1>`__
                  -  `Optional arguments <#optional-arguments_1>`__
                  -  `Example command <#example-command_1>`__
                  -  `Output files <#output-files_1>`__

               -  `Find genes by
                  annotation <#find-genes-by-annotation>`__

                  -  `Required arguments <#required-arguments_2>`__
                  -  `Optional arguments <#optional-arguments_2>`__
                  -  `Example command <#example-command_2>`__
                  -  `Output files <#output-files_2>`__

               -  `Find genes in region <#find-genes-in-region>`__

                  -  `Required arguments <#required-arguments_3>`__
                  -  `Optional arguments <#optional-arguments_3>`__
                  -  `Example input file <#example-input-file>`__
                  -  `Example command <#example-command_3>`__
                  -  `Output files <#output-files_3>`__

            -  `Functional annotations <#functional-annotations>`__

               -  `Show GO <#show-go>`__

                  -  `Required arguments <#required-arguments_4>`__
                  -  `Example commands <#example-commands>`__
                  -  `Output file <#output-file>`__

               -  `Compare GO <#compare-go>`__

                  -  `Required arguments <#required-arguments_5>`__
                  -  `Example command <#example-command_4>`__
                  -  `Output file <#output-file_1>`__

            -  `Homology group
               information <#homology-group-information>`__

               -  `Required arguments <#required-arguments_6>`__
               -  `Optional arguments <#optional-arguments_4>`__
               -  `Example command <#example-command_5>`__
               -  `Output files <#output-files_4>`__

            -  `Sequence alignments <#sequence-alignments>`__
            -  `Matrix files <#matrix-files>`__

               -  `Order matrix <#order-matrix>`__

                  -  `Required argument <#required-argument>`__
                  -  `Optional argument <#optional-argument>`__
                  -  `Example command <#example-command_6>`__
                  -  `Output file <#output-file_2>`__

               -  `Rename matrix <#rename-matrix>`__

                  -  `Required arguments <#required-arguments_7>`__
                  -  `Optional arguments <#optional-arguments_5>`__
                  -  `Example command <#example-command_7>`__
                  -  `Output file <#output-file_3>`__

            -  `Retrieve regions, genomes or
               features <#retrieve-regions-genomes-or-features>`__

               -  `Retrieve regions <#retrieve-regions>`__

                  -  `Required arguments <#required-arguments_8>`__
                  -  `Example command <#example-command_8>`__
                  -  `Example input <#example-input>`__
                  -  `Output file <#output-file_4>`__

               -  `Retrieve features <#retrieve-features>`__

                  -  `Required arguments <#required-arguments_9>`__

               -  `Optional arguments <#optional-arguments_6>`__

                  -  `Example command <#example-command_9>`__
                  -  `Output files <#output-files_5>`__

         -  `Read mapping <../mapping/>`__

         -  `Querying the database <../query/>`__

         -  `Differences between pangenome and
            panproteome <../differences/>`__

         -  `Tutorial <../tutorial/>`__

   .. container:: section wy-nav-content-wrap

      `PanTools manual <..>`__

      .. container:: wy-nav-content

         .. container:: rst-content

            .. container::

               -  `Docs <..>`__ »
               -  Explore the pangenome
               -  

               --------------

            .. container::

               .. container:: section

                  .. rubric:: Explore the pangenome
                     :name: explore-the-pangenome

                  The functionalities on this page allow to actively
                  explore the pangenome.

                  -  Retrieve regions from the pangenome
                  -  Retrieve sequences and functional annotations from
                     homology groups
                  -  Search for genes using a gene name, functional
                     annotation or database node identifier
                  -  Align homology groups or genomic regions
                  -  GO enrichment analysis

                  .. rubric:: Gene locations
                     :name: gene-locations

                  Identify and compare gene clusters of neighbouring
                  genes based on a set of homology groups. First,
                  identifies the genomic position of genes in homology
                  groups, retrieves the order of genes per genome and
                  based on this construct the gene clusters. If homology
                  groups with multiple genomes were selected, the gene
                  cluster composition is compared between genomes. When
                  a ``--phenotype`` is included, gene clusters can be
                  found that only consist of groups of a certain
                  phenotype.

                  For example, 100 groups were predicted as core in a
                  pangenome of 5 genomes. The gene clusters are first
                  identified per genome, whereafter it compares the gene
                  order of one genome to all the other genomes. The
                  result could be 75 groups with genes that are not only
                  homologous but also share their gene neighbourhood.
                  Another example, when accessory (present 2 in to 4
                  genomes) groups are given to this function in
                  combination with a ``--phenotype`` (assigned to only
                  two genomes), the function can return clusters that
                  can only be found in the phenotype members.

                  .. rubric:: Required arguments
                     :name: required-arguments

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--homology-groups``/``-hm`` A text file with
                    homology group node identifiers, seperated by a
                    comma.

                  .. rubric:: Optional arguments
                     :name: optional-arguments

                  | ``--phenotype``/``-ph`` A phenotype name, used to
                    identify gene clusters shared by all phenotype
                    members.
                  | ``--value`` The number of allowed nucleotides
                    between two neighbouring genes (default is 1 MB).
                  | ``--gap-open``/``-go`` When constucting the
                    clusters, allow a number of genes for each cluster
                    that are not originally part of the input groups
                    (default is 0).
                  | ``--core-threshold`` Lower the threshold (%) for a
                    group to be considered core/softcore (default is the
                    total number of genomes found in the groups, not a
                    percentage).
                  | ``--skip``/``-sk`` Exclude a selection of genomes.
                  | ``--reference``/``-ref`` Only include a selection of
                    genomes.
                  | ``--mode ignore-copies`` Duplicated and co-localized
                    genes no longer break up clusters.

                  .. rubric:: Example command
                     :name: example-command

                  ::

                     $ pantools locate_genes -dp tomato_DB -hm phenotype_groups.csv 
                     $ pantools locate_genes -dp tomato_DB -hm unique_groups.csv --value 5000 -go 1 
                     $ pantools locate_genes -dp tomato_DB -hm accessory_groups.csv --core-threshold 95 -go 1

                  .. rubric:: Output files
                     :name: output-files

                  Output files are stored in
                  *database_directory/locate_genes/*

                  -  **gene_clusters_by_position.txt**, the identified
                     gene clusters ordered by their position in the
                     genome.
                  -  **gene_clusters_by_size.txt**, the identified gene
                     clusters ordered from largest to smallest.
                  -  **compare_gene_clusters**, the composition of found
                     gene clusters is compared to the other genomes. For
                     each cluster, it shows which parts match other
                     clusters and which parts do not. The file is not
                     created when homology groups only contain proteins
                     of a single genome (unique).

                  When a ``--phenotype`` is included

                  -  **phenotype_clusters**, homology group node
                     identifiers from phenotype shared and specific
                     clusters.
                  -  **compare_gene_clusters_PHENOTYPE.txt**, the same
                     information as **compare_gene_clusters** but now
                     the gene cluster comparison is only done between
                     phenotype members.

                  --------------

                  .. rubric:: Find genes
                     :name: find-genes

                  .. rubric:: Find genes by name
                     :name: find-genes-by-name

                  Find your genes of interest in the pangenome by using
                  the gene name and extract the nucleotide and protein
                  sequence. To be able to find a gene, every letter of
                  the given input must match a gene name. The search is
                  not case sensitive. Performing a search with 'sonic1'
                  as query will not be able find 'sonic', but is able to
                  find Sonic1, SONIC1 or sOnIc1. Including the
                  ``--mode 1`` argument allows a more relaxed search and
                  using 'sonic' will now also find gene name variations
                  as 'sonic1', 'sonic3' etc..

                  Be aware, for this function to work it is important
                  that genomes are annotated by a method that follows
                  the rules for genetic nomenclature. Gene naming can be
                  inconsistent when different tools are used for genome
                  annotation, making this functionality ineffective.

                  This function is the same as
                  `mlsa_find_genes <../phylogeny/#mlsa>`__ but uses a
                  different output directory. Several warnings (shown in
                  the other manual) can be generated during the seach.
                  These warning are less relevant for this function as
                  the genes are not required to be single
                  copy-orthologous.

                  .. rubric:: Required arguments
                     :name: required-arguments_1

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--name`` One or multiple gene names, seperated by
                    a comma.

                  .. rubric:: Optional arguments
                     :name: optional-arguments_1

                  | ``--skip``/``-sk`` Exclude a selection of genomes.
                  | ``--reference``/``-ref`` Exclude a selection of
                    genomes.
                  | ``--mode extensive`` Perform a more extensive gene
                    search.

                  .. rubric:: Example command
                     :name: example-command_1

                  ::

                     $ pantools find_genes_by_name -dp tomato_DB --name dnaX,gapA,recA
                     $ pantools find_genes_by_name -dp tomato_DB --name gapA --mode extensive

                  .. rubric:: Output files
                     :name: output-files_1

                  Output files are stored in
                  */database_directory/find_genes/by_name/*. For each
                  gene name that was included, a nucleotide and protein
                  and .FASTA file is created with sequences found in all
                  genomes.

                  -  **find_genes_by_name.log**, relevant information
                     about the extracted genes: node identifier, gene
                     location, homology group etc..

                  --------------

                  .. rubric:: Find genes by annotation
                     :name: find-genes-by-annotation

                  Find genes of interest in the pangenome that share a
                  functional annotation node and extract the nucleotide
                  and protein sequence.

                  .. rubric:: Required arguments
                     :name: required-arguments_2

                  ``--database-path``/``-dp`` Path to the database.

                  Requires either **one** of the following arguments

                  | ``--node`` One or multiple identifiers of function
                    nodes (GO, InterPro, PFAM, TIGRFAM), seperated by a
                    comma.
                  | ``--name`` One or multiple function identifiers (GO,
                    InterPro, PFAM, TIGRFAM), seperated by a comma.

                  .. rubric:: Optional arguments
                     :name: optional-arguments_2

                  | ``--skip``/``-sk`` Exclude a selection of genomes.
                  | ``--reference``/``-ref`` Only include a selection of
                    genomes.

                  .. rubric:: Example command
                     :name: example-command_2

                  ::

                     $ pantools find_genes_by_annotation -dp tomato_DB --node 14928,25809
                     $ pantools find_genes_by_annotation -dp tomato_DB --name PF00005,GO:0000160,IPR000683,TIGR02499

                  .. rubric:: Output files
                     :name: output-files_2

                  Output files are stored in
                  */database_directory/find_genes/by_annotation/*. For
                  each function (node) that was included, a nucleotide
                  and protein and .FASTA file is created with sequences
                  from the genes that are connected to the node.

                  -  **find_genes_by_annotation.log**, relevant
                     information about the extracted genes: node
                     identifier, gene location, homology group etc..

                  --------------

                  .. rubric:: Find genes in region
                     :name: find-genes-in-region

                  Find genes of interest in the pangenome that can be
                  (partially) found within a given region (partially).
                  For each found gene, relevant information, the
                  nucleotide sequence and protein sequence is extracted.

                  .. rubric:: Required arguments
                     :name: required-arguments_3

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--regions-file``/``-rf`` A text file containing
                    genome locations with on each line: a genome number,
                    sequence number, begin and end position, separated
                    by a space.

                  .. rubric:: Optional arguments
                     :name: optional-arguments_3

                  ``--mode partial`` Also retrieve genes that only
                  partially overlap the input regions.

                  .. rubric:: Example input file
                     :name: example-input-file

                  Each line must have a genome number, sequence number,
                  begin and end positions that are separated by a space.

                  ::

                     195 1 477722 478426
                     71 10 17346 18056 
                     138 47 159593 160300 

                  .. rubric:: Example command
                     :name: example-command_3

                  ::

                     $ pantools find_genes_in_region -dp tomato_DB -rf regionts.txt
                     $ pantools find_genes_in_region -dp tomato_DB -rf regionts.txt --mode partial 

                  .. rubric:: Output files
                     :name: output-files_3

                  Output files are stored in
                  */database_directory/find_genes/in_region/*. For each
                  region that was included, a nucleotide and protein and
                  .FASTA file is created with sequences from the genes
                  that are found within the region.

                  -  **find_genes_in_region.log**, relevant information
                     about the extracted genes: node identifier, gene
                     location, homology group etc..

                  --------------

                  .. rubric:: Functional annotations
                     :name: functional-annotations

                  The following functions can only be used when any type
                  of functional annotation is `added to the
                  database <../construct/#add-functional-annotations>`__.

                  .. rubric:: Show GO
                     :name: show-go

                  For a selection of '**GO**' nodes, retrieves connected
                  'mRNA' nodes, child and all parent GO terms that are
                  higher in the GO hiercarchy. This function follows the
                  'is_a' relationships of GO each node to their parent
                  GO term until the 'biological process', 'molecular
                  function' or 'cellular location' node is reached. This
                  can be is useful in case InterProScan annotations were
                  included, as these only add the most specifc GO terms
                  of the hierarchy to a sequence.

                  .. rubric:: Required arguments
                     :name: required-arguments_4

                  ``--database-path``/``-dp`` Path to the database

                  Requires either **one** of the following arguments

                  | ``--node`` One or multiple identifiers of 'GO'
                    nodes, seperated by a comma.
                  | ``--name`` One or multiple GO term identifiers,
                    seperated by a comma.

                  .. rubric:: Example commands
                     :name: example-commands

                  ::

                     $ pantools show_go -dp tomato_DB --node 15078,15079
                     $ pantools show_go -dp tomato_DB --name GO:0000001,GO:0000002,GO:0008982

                  .. rubric:: Output file
                     :name: output-file

                  -  **show_go.txt**, information of the selected GO
                     node(s): the connected 'mRNA' nodes, the GO layer
                     below, and all layers above.

                  --------------

                  .. rubric:: Compare GO
                     :name: compare-go

                  Check if and how similar two given GO terms are. For
                  both nodes, follows the 'is_a' relationships up to
                  their parent GO terms until the 'biological process',
                  'molecular function' or 'cellular location' node is
                  reached. After all parent terms are found, the shared
                  GO terms and their location in the hierarchy is
                  reported.

                  .. rubric:: Required arguments
                     :name: required-arguments_5

                  ``--database-path``/``-dp`` Path to the database.

                  Requires either **one** of the following arguments

                  | ``--node`` Two node identifiers of 'GO' nodes,
                    seperated by a comma.
                  | ``--name`` Two GO identifiers, seperated by a comma.

                  .. rubric:: Example command
                     :name: example-command_4

                  ::

                     $ pantools compare_go -dp tomato_DB --name GO:0032775,GO:0006313
                     $ pantools compare_go -dp tomato_DB --node 741487,741488

                  .. rubric:: Output file
                     :name: output-file_1

                  Output files are stored in
                  *database_directory/function/*

                  -  **compare_go.txt**, information of the two GO
                     nodes: the connected 'mRNA' nodes, the GO layer
                     below, all layers above and the shared GO terms
                     between the two nodes.

                  --------------

                  .. rubric:: Homology group information
                     :name: homology-group-information

                  Report all available information of one or multiple
                  homology groups.

                  .. rubric:: Required arguments
                     :name: required-arguments_6

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--homology-groups``/``-hm`` A text file with
                    homology group node identifiers, seperated by a
                    comma

                  .. rubric:: Optional arguments
                     :name: optional-arguments_4

                  | ``--label`` Name of function identifiers from GO,
                    PFAM, InterPro or TIGRAM. To find Phobius (P) or
                    SignalP (S) annotations, include: 'secreted' (P/S),
                    'receptor' (P/S), or 'transmembrane' (P).
                  | ``--name`` One or multiple gene names, seperated by
                    a comma.
                  | ``--skip``/``-sk`` Exclude a selection of genomes.
                  | ``--reference``/``-ref`` Only include a selection of
                    genomes.

                  .. rubric:: Example command
                     :name: example-command_5

                  ::

                     $ pantools group_info -dp yeast_DB -hm core_groups.txt
                     $ pantools group_info -dp yeast_DB -hm core_groups.txt --label GO:0032775,GO:0006313 --name budC,estP

                  .. rubric:: Output files
                     :name: output-files_4

                  Output files are stored in
                  *database_directory/alignments/grouping_v?/groups/*.
                  For each homology group that was included, a
                  nucleotide and protein and .FASTA file is created with
                  sequences found in all genomes.

                  -  **group_info.txt**, relevant information for each
                     homology group: number of copies per genome, gene
                     names, mRNA node identifiers, functions, protein
                     sequence lengths, etc..
                  -  **group_functions.txt**, full description of the
                     functions found in homology groups

                  When function identifiers are included via ``--label``

                  -  **groups_with_function.txt**, homology group node
                     identifiers from groups that match one of the input
                     functions.

                  When gene names are included via ``--name``

                  -  **groups_with_name.txt**, homology group node
                     identifiers from groups that match one of the input
                     gene ames.

                  --------------

                  .. rubric:: Sequence alignments
                     :name: sequence-alignments

                  The manual for PanTools' sequence alignment
                  functionalities moved to a standalone page - `Sequence
                  alignments <../msa/>`__

                  --------------

                  .. rubric:: Matrix files
                     :name: matrix-files

                  Several functions generate tables in a CSV file
                  format. as tables that the following functions can
                  work with. For example, ANI scores, *k*-mer and gene
                  distance used for constructing the Neighbour Joining
                  `phylogenetic trees <../phylogeny/#phylogeny>`__, and
                  the identity and protein sequence similarity tables
                  created by the `alignment
                  functions <./#sequence-alignments>`__.

                  .. rubric:: Order matrix
                     :name: order-matrix

                  Transforms the CSV table to easy to read file by
                  ordering the values in ascending order from low to
                  high or descending order when ``--mode desc`` is
                  included in the command. If phenotype information is
                  included in the header, a seperate file with the range
                  of found values is created for each phenotype. If this
                  information is not present (only genome numbers in the
                  header), use `rename_matrix <./#rename-matrix>`__ to
                  change the headers.

                  .. rubric:: Required argument
                     :name: required-argument

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--input-file``/``-af`` A CSV formatted matrix
                    file.

                  .. rubric:: Optional argument
                     :name: optional-argument

                  | ``--skip``/``-sk`` Skip over the values of a
                    selection of genomes.
                  | ``--reference``/``-ref`` Only include the values
                    from a selection of genomes.
                  | ``--mode asc`` or ``--mode desc`` Order the matrix
                    in ascending or descending order (ascending is
                    default).

                  .. rubric:: Example command
                     :name: example-command_6

                  ::

                     $ pantools order_matrix -dp bacteria_DB -if bacteria_DB/ANI/fastANI/ANI_distance_matrix.csv
                     $ pantools order_matrix -dp bacteria_DB -if bacteria_DB/ANI/fastANI/ANI_distance_matrix.csv --mode desc

                  .. rubric:: Output file
                     :name: output-file_2

                  Output is written to the same directory as the
                  selected input file

                  -  '*old file name*' + '*\_ORDERED*', ordered values
                     of the original matrix file.

                  When phenotype information is present in the header

                  -  '*old file name*' + '*\_PHENOTYPE*', range of
                     values per phenotype.

                  --------------

                  .. rubric:: Rename matrix
                     :name: rename-matrix

                  Rename the headers (first row and leftmost column) of
                  CSV formatted matrix files. If no ``--phenotype`` is
                  included, headers are changed to only contain genome
                  numbers.

                  .. rubric:: Required arguments
                     :name: required-arguments_7

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--input-file``/``-af`` a matrix file with
                    numerical values.

                  .. rubric:: Optional arguments
                     :name: optional-arguments_5

                  | ``--phenotype``/``-ph`` A phenotype name, used to
                    include phenotype information into the headers.
                  | ``--skip``/``sk`` Exclude a selection of genomes
                    from the new matrix file. ``--reference``/``-ref``
                    Only include a selection of genomes in the new
                    matrix file.
                  | ``--mode no-numbers`` Exclude genome numbers from
                    the headers.

                  .. rubric:: Example command
                     :name: example-command_7

                  ::

                     $ pantools rename_matrix -dp pecto_DB -phenotype species -if pecto_DB/ANI/fastANI/ANI_distance_matrix.csv

                  .. rubric:: Output file
                     :name: output-file_3

                  Output is written to the same directory as the
                  selected input file.

                  -  '*old file name*' + '*\_RENAMED*', the original
                     matrix file with changed headers.

                  --------------

                  .. rubric:: Retrieve regions, genomes or features
                     :name: retrieve-regions-genomes-or-features

                  The two following functions allow users to retrieve
                  genomic regions from the pangenome.

                  .. rubric:: Retrieve regions
                     :name: retrieve-regions

                  Retrieve the full genome sequence or genomic regions
                  from the pangenome.

                  .. rubric:: Required arguments
                     :name: required-arguments_8

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--regions-file``/``-rf`` A text file containing
                    genome locations with on each line: a genome number,
                    sequence number, begin and end positions separated
                    by a space.

                  .. rubric:: Example command
                     :name: example-command_8

                  ::

                     $ pantools retrieve_regions -dp pecto_DB -rf regions.txt

                  .. rubric:: Example input
                     :name: example-input

                  To extract:

                  -  Complete genome - Include a genome number
                  -  An entire sequence - Include a genome number with
                     sequence number
                  -  A genomic region - Include a genome number,
                     sequence number, begin and end positions that are
                     separated by a space. Place a minus symbol behind
                     the regions to extract the reverse complement
                     sequence of the region.

                  ::

                     1
                     1 1 
                     1 1 1 10000
                     1 1 1000 1500 -
                     195 1 477722 478426
                     71 10 17346 18056 -
                     138 47 159593 160300 -

                  .. rubric:: Output file
                     :name: output-file_4

                  A single FASTA file is created for all given locations
                  and is stored in the database directory.

                  .. rubric:: Retrieve features
                     :name: retrieve-features

                  To retrieve the sequence of annotated features from
                  the pangenome.

                  .. rubric:: Required arguments
                     :name: required-arguments_9

                  | ``--database-path``/``-dp`` Path to the database.
                  | ``--feature-type`` or ``-ft`` The feature name; for
                    example 'gene', 'mRNA', 'exon', 'tRNA', etc.

                  .. rubric:: Optional arguments
                     :name: optional-arguments_6

                  Use one of the following arguments to limit the
                  seqeunce retrieval to a selection of genomes.

                  | ``--skip``/``-sk`` Exclude a selection of genomes.
                  | ``--reference``/``-ref`` Only include a selection of
                    genomes.

                  .. rubric:: Example command
                     :name: example-command_9

                  ::

                     $ pantools retrieve_features -dp pecto_DB --feature-type gene
                     $ pantools retrieve_features -dp pecto_DB --ft mRNA

                  .. rubric:: Output files
                     :name: output-files_5

                  For each genome a FASTA file containing the retrieved
                  features will be stored in the database directory. For
                  example, genes.1.fasta contains all the genes
                  annotated in genome 1.

                  --------------

            .. container:: rst-footer-buttons

               `Next <../mapping/>`__ `Previous <../msa/>`__

            --------------

            .. container::

.. container:: rst-versions

   `« Previous <../msa/>`__ `Next » <../mapping/>`__
