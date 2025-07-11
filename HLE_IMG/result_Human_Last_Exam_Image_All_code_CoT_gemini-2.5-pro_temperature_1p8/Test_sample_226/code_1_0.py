import collections

# Based on the visual analysis of the immunohistochemistry images, there is a clear
# decrease in the number of APT1-positive cells (brown stain) in the brains of
# patients with Parkinson's Disease (PD) and Parkinson's Disease with dementia (PDD)
# compared to healthy controls.

# Let's represent the quantitative data from answer choice A to verify this trend.
# The data is given as mean ± standard deviation in cells per mm^2.

# Data from statement A:
cell_counts = {
    "control": {"mean": 679.6, "std_dev": 59.32},
    "PD": {"mean": 302.1, "std_dev": 111.5},
    "PDD": {"mean": 283.2, "std_dev": 42.26}
}

print("Quantitative analysis of APT1 immunopositive cells based on statement A:")
for group, data in cell_counts.items():
    print(f"- {group.upper()} group: {data['mean']} ± {data['std_dev']} cells/mm²")

print("\nConclusion based on this data:")
# Compare control to the disease groups
control_mean = cell_counts['control']['mean']
pd_mean = cell_counts['PD']['mean']
pdd_mean = cell_counts['PDD']['mean']

if control_mean > pd_mean and control_mean > pdd_mean:
    print(f"The number of cells in the control group ({control_mean}) is substantially higher than in the PD ({pd_mean}) and PDD ({pdd_mean}) groups.")

# Compare PD and PDD groups
if abs(pd_mean - pdd_mean) < (pd_mean * 0.1): # check if difference is less than 10% of the PD value
     print(f"The cell counts for PD ({pd_mean}) and PDD ({pdd_mean}) are relatively similar to each other.")

print("\nThis quantitative relationship (Control > PD ≈ PDD) aligns with the visual evidence presented in the images.")
print("Therefore, the statement providing this data is the most plausible.")
