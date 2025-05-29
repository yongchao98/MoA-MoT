# Define the truth values for each person
andree_truth = False  # Andree lies
amberly_truth = andree_truth  # Amberly says Andree tells the truth
gwenn_truth = amberly_truth  # Gwenn says Amberly tells the truth
willian_truth = not gwenn_truth  # Willian says Gwenn lies
alexis_truth = willian_truth  # Alexis says Willian tells the truth

# Output whether Alexis tells the truth
print(alexis_truth)