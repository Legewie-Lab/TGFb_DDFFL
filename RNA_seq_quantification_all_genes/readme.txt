The CSV file contains the quantification of mRNA expression (RPKM) for all genes utilized in this study. Below are details about the columns:
Column Descriptions

    Gene Identifiers:
    The rows typically contain gene names or IDs.

    RPKM Values:
    RNA-seq data is quantified using RPKM (Reads Per Kilobase of transcript per Million mapped reads). This normalization accounts for gene length and sequencing depth.

    Columns R1, R2, R3:
    Indicate biological replicates.
        R1, R2, R3: The replicate number.
        ctl: Indicates control samples.
        low/high: Indicates TGFβ dosage (low or high concentration).

Example Column Naming Convention:

    R1_ctl: Replicate 1, control condition.
    R2_high: Replicate 2, high TGFβ dosage.
    R3_low: Replicate 3, low TGFβ dosage
