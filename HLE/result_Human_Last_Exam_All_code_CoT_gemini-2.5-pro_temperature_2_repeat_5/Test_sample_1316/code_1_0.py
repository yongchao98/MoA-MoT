# A program to calculate the chromosome number of an exceptional X0 Drosophila male.
# In Drosophila, there are 3 pairs of autosomes.
num_autosome_pairs = 3

# Calculate the total number of autosomes.
num_autosomes = num_autosome_pairs * 2

# The exceptional male is X0, meaning he has one X chromosome and no second sex chromosome (Y).
num_sex_chromosomes = 1

# The total number of chromosomes is the sum of autosomes and sex chromosomes.
total_chromosomes = num_autosomes + num_sex_chromosomes

print("Analysis of the Exceptional X0 Male's Chromosomes:")
print(f"Number of Autosomes: {num_autosomes}")
print(f"Number of Sex Chromosomes (X0): {num_sex_chromosomes}")
print("\nFinal equation for the total chromosome count:")
print(f"{num_autosomes} + {num_sex_chromosomes} = {total_chromosomes}")
print(f"\nThis is an aneuploid condition with a total of {total_chromosomes} chromosomes, compared to the normal 8.")