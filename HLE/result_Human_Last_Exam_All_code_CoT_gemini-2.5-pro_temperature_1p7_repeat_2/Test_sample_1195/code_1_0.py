# Define the components for the autosomal cross (su-v+/su-v x su-v+/su-v)
# Genotypes and their counts in a standard monohybrid cross out of 4 total offspring.
homozygous_dominant = 1  # su-v+/su-v+
heterozygous = 2         # su-v+/su-v
homozygous_recessive = 1 # su-v/su-v

total_offspring = homozygous_dominant + heterozygous + homozygous_recessive

# Determine the phenotype based on the su-v gene, as all F2 offspring have the
# genetic background for vermilion eyes (X(v)X(v) or X(v)Y).
# The recessive su-v/su-v genotype suppresses the vermilion phenotype, resulting in wild-type eyes.
# The su-v+ allele is dominant and does not suppress.

# Offspring that are NOT suppressed (show vermilion phenotype)
vermilion_count = homozygous_dominant + heterozygous

# Offspring that ARE suppressed (show wild-type phenotype)
wild_type_count = homozygous_recessive

# Print the step-by-step logic
print("Step 1: All F1 offspring are heterozygous for the suppressor gene (su-v+/su-v) and have the vermilion genotype (X(v)X(v) or X(v)Y).")
print("Step 2: The F1 cross is X(v)X(v); su-v+/su-v  x  X(v)Y; su-v+/su-v.")
print("Step 3: Since all F2 offspring have the vermilion genotype, the phenotype is determined by the su-v gene cross.")
print("Step 4: The cross su-v+/su-v x su-v+/su-v yields the following ratio of genotypes:")
print(f"  - {homozygous_dominant}/{total_offspring} are su-v+/su-v+ (Non-suppressed -> Vermilion)")
print(f"  - {heterozygous}/{total_offspring} are su-v+/su-v (Non-suppressed -> Vermilion)")
print(f"  - {homozygous_recessive}/{total_offspring} are su-v/su-v (Suppressed -> Wild-type)")
print("\nStep 5: Calculating the final phenotypic ratio...")
print(f"The number of vermilion offspring is {homozygous_dominant} + {heterozygous} = {vermilion_count}")
print(f"The number of wild-type offspring is {wild_type_count}")

print("\n-------------------------------------------")
print(f"Final F2 Phenotypic Ratio:")
print(f"{vermilion_count}/{total_offspring} vermilion : {wild_type_count}/{total_offspring} wild-type")
print("-------------------------------------------")