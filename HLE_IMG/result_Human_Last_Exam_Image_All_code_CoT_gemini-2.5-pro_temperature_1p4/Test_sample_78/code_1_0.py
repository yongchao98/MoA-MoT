# Step 1: Define the counts for each pattern observed on the gel.
# Pattern 1: Homozygous type 1 (e.g., single high band)
# Pattern 2: Heterozygous (e.g., two bands)
# Pattern 3: Homozygous type 2 (e.g., single low band)
count_pattern1 = 2
count_pattern2 = 8
count_pattern3 = 6

# Step 2: Assign patterns to genotypes based on genetic principles.
# The pattern with two bands (or a mix) is the heterozygote.
# The observed count for this pattern matches the Mendelian expectation (1/2 * 16 = 8).
num_heterozygotes = count_pattern2

# The other two patterns are the homozygotes. Their counts are 2 and 6.
homozygote_counts = [count_pattern1, count_pattern3]

# Step 3: Use a biological assumption to distinguish between wild-type and mutant homozygotes.
# Assumption: The homozygous mutant genotype is likely to have reduced viability,
# resulting in a lower count than expected. We assign the smaller count to the mutants.
num_homozygous_mutants = min(homozygote_counts)
num_homozygous_wild_type = max(homozygote_counts)

# Step 4: Finalize the counts for A, B, and C.
# A = # of homozygous wild-type larvae
# B = # of heterozygotes
# C = # of homozygous mutants
A = num_homozygous_wild_type
B = num_heterozygotes
C = num_homozygous_mutants

# Print the explanation and the final answer.
print("Based on the analysis of the gel and principles of Mendelian genetics:")
print(f"Number of homozygous wild-type larvae (A) = {A}")
print(f"Number of heterozygous larvae (B) = {B}")
print(f"Number of homozygous mutant larvae (C) = {C}")
print("\nThe final answer in the format A/B/C is:")
print(f"{A}/{B}/{C}")