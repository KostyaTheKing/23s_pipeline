from pathlib import Path
import pandas as pd

# Structure of the resulting table.
resulting_dataframe = pd.DataFrame(
    columns=[
        "genome_name",  # A name of the genome
        "2611",  # A letter present at site 2611
        "2611_num",  # Number of copies
        "2611_reads_%",  # % of reads thath support SNP
        "2058",
        "2058_num",
        "2058_reads_%",
        "2059",
        "2059_num",
        "2059_reads_%",
    ]
)


def name_replacer(arr, file_path, index):
    """
    A function to replace genome name in vcf file for genome name in fasta.
    """
    arr[index] = file_path.stem
    return arr


def coord_replacer(arr, coords, index):
    """
    A function to replace big coordinates to more readable E. coli coordinates.
    """
    for old_coord, standard_coord in coords.items():
        if old_coord in arr[index]:
            arr[index] = arr[index].replace(old_coord, standard_coord)
    return arr


# WHO_F_2024 SNP coordinates to E. coli coordinates
ecoli_coordinates = {"1063174": "2058", "1063175": "2059", "1063727": "2611"}


def vcf_iterator(file_list, results, ecoli_coordinates):
    """
    A function to iterate through vcf files and write the results in a table.
    """
    for i in file_list:
        # Reading vcf files and looking for the relevant lines.
        with open(i, "r") as f:
            relevant_lines = [
                line
                for line in f
                if sum(j in line for j in ecoli_coordinates.keys())
                > 0  # If the relevant coordinate is present in a line, then save it in a list.
            ]
        if (
            len(relevant_lines) == 0
        ):  # If nothing is found, than there are no polymorphisms, i.e. write default values.
            # resulting_dataframe.loc[len(resulting_dataframe),:] = [i.stem + '.fas', 0, 0, 0, 0, 0, 0]
            resulting_dataframe.loc[len(resulting_dataframe), :] = [
                i.stem,
                "C",
                4,
                100,
                "A",
                4,
                100,
                "A",
                4,
                100,
            ]
        else:
            # Create a line with default values and then edit it later.
            res_of_analysis = [i.stem, "C", 4, 100, "A", 4, 100, "A", 4, 100]

            # Split read lines by tab.
            relevant_lines = map(lambda x: x.split("\t"), relevant_lines)

            # Replace genome name to the user provided name
            relevant_lines = map(lambda x: name_replacer(x, i, 0), relevant_lines)

            # Extract important information from relevant lines
            relevant_lines = map(
                lambda x: [
                    x[0],  # Genome name
                    x[1],  # Coordinate
                    x[3],  # Default letter at SNP site
                    x[4],  # Found letter at SNP site
                    "",  # Empty place for genotype
                    x[9]
                    .split(":")[2]
                    .split(
                        ","
                    ),  # Total amount of reads for each allele comma separated.
                ],
                relevant_lines,
            )

            # Replace amount of reads with % of reads that support SNP
            relevant_lines = map(
                lambda x: x[:-1]
                + [str(float(x[-1][1]) / (float(x[-1][0]) + float(x[-1][1])) * 100)],
                relevant_lines,
            )

            # Round % of reads to the nearest quartile to determine SNP copy number
            quartile_list = [0, 25, 50, 75, 100]
            round_to_nearest_quartile = lambda number, quartile_list: min(
                quartile_list, key=lambda y: abs(y - float(number))
            )
            relevant_lines = map(
                lambda x: x[:4]
                + [
                    str(
                        quartile_list.index(
                            round_to_nearest_quartile(x[-1], quartile_list)
                        )
                    )
                ]
                + x[5:],
                relevant_lines,
            )

            # Replace coordinates with E. coli ones.
            relevant_lines = map(
                lambda x: coord_replacer(x, ecoli_coordinates, 1), relevant_lines
            )

            relevant_lines = [i for i in relevant_lines if i[4] != "0"]

            if not relevant_lines:
                resulting_dataframe.loc[len(resulting_dataframe), :] = res_of_analysis
                continue

            # Write data to the DataFrame
            for found_snps in relevant_lines:
                if "2058" in found_snps:
                    res_of_analysis[4] = found_snps[3]
                    res_of_analysis[5] = int(found_snps[4])
                    res_of_analysis[6] = float(found_snps[-1])
                if "2059" in found_snps:
                    res_of_analysis[7] = found_snps[3]
                    res_of_analysis[8] = int(found_snps[4])
                    res_of_analysis[9] = float(found_snps[-1])
                if "2611" in found_snps:
                    res_of_analysis[1] = found_snps[3]
                    res_of_analysis[2] = int(found_snps[4])
                    res_of_analysis[3] = float(found_snps[-1])
                resulting_dataframe.loc[len(resulting_dataframe), :] = res_of_analysis

    resulting_dataframe.to_excel(results, index=False)


if __name__ == "__main__":
    # Path to the folder with vcf files
    file_list_path = Path("./rrn_genome_variants")

    # Read file names in a folder
    file_list = file_list_path.iterdir()

    # Sort files in numerical or alphabetical order
    try:
        file_list = sorted(file_list, key=lambda x: int(x.stem))
    except:
        file_list = sorted(file_list, key=lambda x: x.stem)

    # Путь к файлу с результатами анализа vcf файлов
    results = Path("./results.xlsx")

    # Remove an old file if present
    results.unlink(missing_ok=True)

    # Launch vcf_iterator
    vcf_iterator(
        file_list=file_list, results=results, ecoli_coordinates=ecoli_coordinates
    )