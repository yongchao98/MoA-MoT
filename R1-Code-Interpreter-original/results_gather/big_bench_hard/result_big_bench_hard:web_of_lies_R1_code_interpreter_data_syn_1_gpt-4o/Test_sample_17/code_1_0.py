# Define the truth values for each person
ka_truth = False  # Ka lies
andree_truth = not ka_truth  # Andree says Ka tells the truth
audrie_truth = not andree_truth  # Audrie says Andree lies
antwan_truth = audrie_truth  # Antwan says Audrie tells the truth
millie_truth = antwan_truth  # Millie says Antwan tells the truth

# Print the result for Millie's truthfulness
print(millie_truth)