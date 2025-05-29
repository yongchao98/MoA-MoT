# Define the truth values
amberly_truth = False
shenna_truth = amberly_truth  # Shenna says Amberly tells the truth
willian_truth = shenna_truth  # Willian says Shenna tells the truth
raymond_truth = not willian_truth  # Raymond says Willian lies
sal_truth = not raymond_truth  # Sal says Raymond lies

# Output whether Sal tells the truth
print(sal_truth)