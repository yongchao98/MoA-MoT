# This script explains the genetic logic behind the formation of the
# exceptional X0 male in Bridges' Drosophila experiments.

# 1. Define the parental cross and the observed exceptional offspring.
#    Let's use '+' for red eyes/wild-type and 'w' for white eyes.
#    The male also carries the 'm' (miniature wing) allele.
parent_female = "XwXw"  # White-eyed female
parent_male = "X+mY"   # Red-eyed, miniature-winged male
exceptional_offspring = "X+m0" # The observed male offspring

print(f"Context: A cross between a {parent_female} female and a {parent_male} male.")
print(f"Unexpected Offspring Found: A male with genotype {exceptional_offspring}.")
print("-" * 60)

# 2. Deconstruct the offspring's genotype to find the parental gametes.
print("To form an X+m0 zygote, we need two gametes:")
gamete_from_father = "X+m"
gamete_from_mother = "0"
print(f" - A sperm from the father carrying the '{gamete_from_father}' chromosome.")
print(f" - An egg from the mother carrying no sex chromosome (a '{gamete_from_mother}' gamete).")
print("-" * 60)

# 3. Explain how the abnormal '0' egg is formed.
print("A '0' egg is the result of non-disjunction in the female during meiosis.")
print(f"The mother starts with the genotype {parent_female}.")
print("During Meiosis I, the two Xw homologous chromosomes should be pulled to opposite poles.")
print("In non-disjunction of Meiosis I, both Xw chromosomes go to one pole, and none go to the other.")
print("This results in one egg precursor with 'XwXw' and another that becomes the '0' egg.")
print("-" * 60)

# 4. Display the final equation of fertilization for this specific offspring.
print("Final Equation:")
# The print statement below shows each part of the "equation" as requested.
print(f"Mother's abnormal egg: {gamete_from_mother}")
print(f"Father's normal sperm: {gamete_from_father}")
print(f"Resulting Zygote: {exceptional_offspring}")
print(f"\nEquation: '{gamete_from_mother}' (from Female) + '{gamete_from_father}' (from Male) -> '{exceptional_offspring}' (Offspring)")
print("-" * 60)
print("This event, non-disjunction of the X chromosome in female meiosis I, explains the presence of the unexpected X0 male.")
