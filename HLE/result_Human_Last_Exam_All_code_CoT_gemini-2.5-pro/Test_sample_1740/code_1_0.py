# In humans, the let-7 family of microRNAs is encoded by several genes.
# We will sum the number of genes for each member of the family.

# Dictionary mapping let-7 family members to the number of genes encoding them.
let7_genes = {
    "let-7a": 3,
    "let-7b": 1,
    "let-7c": 1,
    "let-7d": 1,
    "let-7e": 1,
    "let-7f": 2,
    "let-7g": 1,
    "let-7i": 1,
    "miR-98": 1, # miR-98 is a member of the let-7 family
}

# Get the counts for each member
gene_counts = list(let7_genes.values())

# Calculate the total number of genes
total_members = sum(gene_counts)

# Create the equation string for printing
# The string.join() method requires strings, so we convert numbers to strings
equation_parts = [str(count) for count in gene_counts]
equation_str = " + ".join(equation_parts)

# Print the result
print(f"The total number of human let-7 family members is the sum of genes for each type:")
print(f"{equation_str} = {total_members}")