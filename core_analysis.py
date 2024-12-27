# Packages and data
import os
from copy import deepcopy
import polars as pl
import time
import statistics

# # # # # LOG # # # #
# 1. For some reason, pl.concat does not accept 2 lazyframe otherwise memory leakage → program does not reach its
# optimum
# 2. Package pl is not stable enough → amenable to future adjustment
# 3. Change the directory to cater with all computers

# Directory and files

os.chdir("/Users/xuzifeng/Desktop/bit150/proj2")


def allele_tab_2_dataframe(filename, snp_id) -> pl.lazyframe:
    os.chdir("/Users/xuzifeng/Desktop/bit150/proj2/dbsnps")

    with open(filename, 'r') as file:
        lines = file.readlines()[12:]
        file.close()
    with open('../intermediate_files/allele_intermediate.tsv', 'w') as intermediate:
        intermediate.writelines(lines)

    df = (pl.scan_csv(source="/Users/xuzifeng/Desktop/bit150/proj2/intermediate_files/allele_intermediate.tsv",
                      separator='\t')
          .drop(['BioProject ID', 'BioSample ID', 'Group'])
          .with_columns(pl.lit(snp_id).alias('snp_id')))
    # .explain(optimized=True)
    # .collect())

    os.chdir("/Users/xuzifeng/Desktop/bit150/proj2")

    return df


def one_single_nt_in_single_alt_allele_cell(frame: pl.LazyFrame) -> bool:
    single_element = True

    alt_alleles = frame.select(pl.col('Alt Allele')).unique().collect().to_series()
    for alt in alt_alleles:
        if len(alt.split(',')) != 1:
            single_element = False
            return single_element

    return single_element


def table_normalizer(dont_need_to_be_normalized: bool, target_table: pl.LazyFrame) -> pl.LazyFrame:
    def reshape_and_cast(normal: pl.LazyFrame) -> pl.LazyFrame:
        normal = (normal.with_columns(pl.col("Ref Allele").str.split("=").list.get(0).alias("Ref"),
                                      pl.col("Ref Allele").str.split("=").list.get(1).alias("Freq(Ref)"),
                                      pl.col("Alt Allele").str.split("=").list.get(0).alias("Alt"),
                                      pl.col("Alt Allele").str.split("=").list.get(1).alias("Freq(Alt)"))
                  .drop(["Alt Allele", 'Ref Allele'])
                  .cast({"Freq(Ref)": pl.Float64, "Freq(Alt)": pl.Float64, "Samplesize": pl.Float64})
                  )
        return normal

    if dont_need_to_be_normalized:
        return reshape_and_cast(target_table)

    else:
        # Initialize by mapping. Column name : []
        normalized_table = pl.LazyFrame(data={col: [] for col in target_table.collect_schema().names()})

        for num_rows in range(len(target_table.collect())):
            curr_row = deepcopy(target_table.limit(num_rows).last())
            is_good = one_single_nt_in_single_alt_allele_cell(curr_row)

            if is_good:
                normalized_table = pl.concat([curr_row, normalized_table])

            else:
                alt_alleles = curr_row.select(pl.col("Alt Allele")).collect()[0, 0].split(", ")

                for alle in alt_alleles:
                    # Create a copy for each allele, and append it to the lazyframe
                    lf = pl.LazyFrame({"Alt Allele": alle})
                    curr_row = curr_row.update(lf, how='full')
                    normalized_table = pl.concat([curr_row, normalized_table])

        return reshape_and_cast(normalized_table)


def allele_analysis(normalized_table: pl.LazyFrame) -> pl.lazyframe:
    return (normalized_table.group_by(["Population", "Ref", "Alt", 'snp_id'], maintain_order=True).
            agg([pl.col("Samplesize"), pl.col("Freq(Alt)"), pl.col("Freq(Ref)"), pl.col("#Study")])
            .with_columns(
        ((pl.col("Samplesize") * pl.col("Freq(Alt)")).list.sum() / pl.col("Samplesize").list.sum()).alias(
            'Weighted_Freq(Alt)'))
            .with_columns(
        ((pl.col("Samplesize") * pl.col("Freq(Ref)")).list.sum() / pl.col("Samplesize").list.sum()).alias(
            'Weighted_Freq(Ref)'))
            .with_columns(pl.col("Samplesize").list.sum().alias("sample_size"))
            .drop(['Freq(Alt)', 'Freq(Ref)', 'Samplesize']))


def get_merged_table():
    # # # Get Merged Allele Table
    os.chdir("/Users/xuzifeng/Desktop/bit150/proj2/dbsnps")  # move to large merge function later

    all_snp_files = os.listdir()
    all_snp_files = [snp_f for snp_f in all_snp_files if snp_f.startswith('rs')]
    primer_file = all_snp_files[0]
    primer_snp = primer_file.split('_')[0]
    pre_primer = allele_tab_2_dataframe(primer_file, primer_snp)
    primer = allele_analysis(table_normalizer(one_single_nt_in_single_alt_allele_cell(pre_primer), pre_primer))

    output = primer.collect()
    passed_files = []

    try:
        for file in all_snp_files[1:]:
            snp = file.split('_')[0]
            passed_files.append(file)

            # For some reason, both of the files have to be dataframe otherwise memory leak

            target = allele_tab_2_dataframe(file, snp)
            n = allele_analysis(
                table_normalizer(
                    one_single_nt_in_single_alt_allele_cell(target), target)).collect()

            output = pl.concat([output, n])
    except Exception as e:
        print(passed_files.pop())
        print(len(passed_files))
        print(e)

    # # # Get merged table: join with the gene table

    os.chdir("/Users/xuzifeng/Desktop/bit150/proj2")

    return output


# Time
t = []

for i in range(10):
    start = time.time()
    get_merged_table()
    end = time.time()
    t.append(end - start)

print(statistics.mean(t))
