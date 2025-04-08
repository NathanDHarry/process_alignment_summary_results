import pandas as pd
import os
import re
import warnings
import seaborn as sns
import matplotlib.pyplot as plt

# Define the folder containing the summary files created by samtools 'flagstat'.
path_to_files = r"C:\Users\ndhar\Downloads\popgen_alignment_summaries"

# Create a list to store the summaries as DataFrames (from pandas).
summary_data = []

# Loop through all the summary files in the folder.
for file_name in os.listdir(path_to_files):
    if file_name.endswith(".tsv"):

        # Print each .tsv file's name.
        print(f"File: {file_name}")

        # Get the file's global path.
        file_path = os.path.join(path_to_files, file_name)

        # Extract the sample name from the name of the file.
        sample_name = re.search(r"^(.*?)_L001", file_name)
        sample_name = sample_name.group(1) if sample_name else None

        # Split the file name into its component parts.
        name_components = file_name.split("_")

        # Extract the reference genome for the alignment that produced the file.
        reference = None
        for component in name_components:
            if component in ["bay", "lb", "sg"]:
                reference = component
                break
        if reference is None:
            warnings.warn(f"No match found for 'bay', 'lb', or 'sg' in file: {file_name}")

        # Read the contents of the file into a DataFrame object (using Pandas).
        data = pd.read_csv(file_path, sep="\t", header=None)

        # Extract specific values from each DataFrame and save them in variables.
        try:
            total_reads = int(data[data[2] == "total (QC-passed reads + QC-failed reads)"].iloc[0, 0])
            mapped_reads = int(data[data[2] == "mapped"].iloc[0, 0])
            mapped_percent = round((mapped_reads / total_reads) * 100, ndigits=2)

        except IndexError:
            total_reads = None
            mapped_reads = None
            mapped_percent = None

        # Create a dictionary for each file's extracted attributes.
        data = {
            "FileName": file_name,
            "SampleName": sample_name,
            "Reference": reference,
            "TotalReads": total_reads,
            "MappedReads": mapped_reads,
            "MappedPercent": mapped_percent
        }

        # Add the dictionary to the list if all the attributes are present.
        if all(value is not None for value in data.values()):
            summary_data.append(data)
        else:
            warnings.warn(f"File '{file_name}' has missing values: {data}")

# Create a single DataFrame from the list of dictionaries
summary_df = pd.DataFrame(summary_data)


# Generate a boxplot of 'mapped_percent' for each category of 'reference' (bay, lb, and sg).
plt.figure(figsize=(10, 6))

# Create a custom palette
custom_palette = {"bay": "forestgreen", "lb": "darkorange", "sg": "lightseagreen"}

sns.boxplot(data=summary_df,
            x="Reference",
            y="MappedPercent",
            order=["bay", "lb", "sg"],
            hue="Reference",
            palette=custom_palette,
            dodge=False,
            fliersize=0)

sns.stripplot(
    data=summary_df,
    x="Reference",
    y="MappedPercent",
    order=["bay", "lb", "sg"],
    hue="Reference",
    palette=custom_palette,
    dodge=False,  # Ensures jitter is applied to the same location as the boxplot
    jitter=0.3,  # Adds jitter to spread out the points horizontally
    alpha=0.7,    # Adjust transparency
    linewidth=0.5,  # Outline for better visibility
    edgecolor="auto"  # Define the edge color for points
)

# Add labels and title to the boxplot
plt.title("Distribution of pop-gen data mapping rates for three genomes.", fontsize=14)
plt.xlabel("Reference Genome", fontsize=12)
plt.ylabel("Mapped Percent (%)", fontsize=12)

# Make the Reference labels (X-axis capitalised)
plt.xticks(ticks=plt.xticks()[0], labels=["BAY", "LB", "SG"], fontsize=10)

# Adjust legend to avoid duplication
plt.legend([], [], frameon=False)  # Hide the stripplot legend if hue is repetitive

# Save the plot in the directory containin the files.
output_path = r"C:\Users\ndhar\Downloads\popgen_reads_mapped_percent_boxplot.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Define a function to identify outlier samples (based on mapping rate) using the IQR method.
def find_outliers_by_reference(df, reference_column, value_column, sample_column):
    outliers_dict = {}
    for ref in df[reference_column].unique():

        # Filter the DataFrame to only those entries for one reference at a time.
        reference_data = df[df[reference_column] == ref]

        # Calculate Q1 (25th percentile)
        q1 = reference_data[value_column].quantile(0.25)
        q3 = reference_data[value_column].quantile(0.75)
        iqr = q3 - q1  # Inter-quartile range

        # Determine the lower and upper bounds for outliers
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr

        # Find outliers
        ref_outliers = reference_data[(reference_data[value_column] < lower_bound) | (reference_data[value_column] > upper_bound)]

        outlier_names = ref_outliers.assign(OutlierName=ref_outliers[value_column].astype(str) + "% " + ref_outliers[sample_column])
        outlier_names = outlier_names.set_index(sample_column)["OutlierName"].to_dict()

        # Store sample names of outliers in the dictionary
        outliers_dict[ref] = outlier_names

    return outliers_dict

# Use the function to identify outliers for 'mapped_percent'.
outliers = find_outliers_by_reference(summary_df, reference_column="Reference", value_column="MappedPercent", sample_column="SampleName")

# Write outliers to a tab-separated file
output_path = r"C:\Users\ndhar\Downloads\popgen_reads_mapped_percent_outliers_by_reference.tsv"
with open(output_path, "w") as f:
    # Write header
    f.write("\t".join(outliers.keys()) + "\n")

    # Write rows of outliers
    outlier_lists = {ref: list(samples.values()) for ref, samples in outliers.items()}  # Convert values to lists
    max_outliers = max(len(outliers_list) for outliers_list in outlier_lists.values())  # Get the maximum number of outliers

    for i in range(max_outliers):
        row = []
        for reference in outlier_lists.keys():
            # Add the outlier if it exists, otherwise add an empty string
            row.append(outlier_lists[reference][i] if i < len(outlier_lists[reference]) else "")
        f.write("\t".join(row) + "\n")

