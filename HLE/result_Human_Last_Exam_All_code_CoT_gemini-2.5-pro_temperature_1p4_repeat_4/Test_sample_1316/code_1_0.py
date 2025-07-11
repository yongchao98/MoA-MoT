# Misión Bridges' Drosophila Experiment Analysis

# 1. Define the Parental Genotypes (based on the classic experiment)
# Let W = red eyes (dominant) and w = white eyes (recessive). The gene is on the X chromosome.
# Female parent is white-eyed.
female_genotype = "XwXw"
# Male parent is red-eyed.
male_genotype = "XWY"

print("Parental Cross:")
print(f"Female (white-eyed): {female_genotype}")
print(f"Male (red-eyed): {male_genotype}\n")

# 2. Normal Meiosis and Offspring
# Female produces only Xw gametes.
# Male produces XW and Y gametes.
normal_female_offspring = "XWXw" # Red-eyed female
normal_male_offspring = "XwY" # White-eyed male

print("Expected Normal Offspring:")
print(f"Daughters: {normal_female_offspring} (Red-eyed)")
print(f"Sons: {normal_male_offspring} (White-eyed)\n")

# 3. Analyze the Unexpected Offspring
# The problem states an unexpected male with red eyes and an X0 karyotype.
# (The miniature wing trait is also X-linked and would follow the same inheritance pattern as the red-eye allele from the father).
unexpected_male_karyotype = "X0"
unexpected_male_phenotype = "Red eyes"

# 4. Trace the Chromosomes of the Unexpected Male
# To be male with an X0 karyotype, the fly received one X chromosome and no Y chromosome.
# To have red eyes, this single X chromosome must carry the 'W' allele.
# The father (XWY) is the only parent with the 'W' allele.
# Therefore, the sperm from the father must have been 'XW'.
sperm_contribution = "XW"
# The mother must have contributed no sex chromosome.
egg_contribution = "0 (nullo-X)"

print("Analysis of Unexpected X0 Red-Eyed Male:")
print(f"This male's genotype must be {sperm_contribution}{egg_contribution}.")
print(f"It inherited the {sperm_contribution} chromosome from the father.")
print(f"It inherited a gamete with no sex chromosome ({egg_contribution}) from the mother.\n")

# 5. Determine the Cause in the Mother (Female Gametogenesis)
# A '0' (nullo-X) egg results from a non-disjunction event during meiosis in the female.
# Let's analyze the two possibilities for non-disjunction in the XwXw female.

print("Mechanism of Non-disjunction in the Female (XwXw):")
# Option A: Non-disjunction in Meiosis I
print("\nScenario A: Non-disjunction in Meiosis I")
print("Step 1: The homologous pair of Xw chromosomes fails to separate.")
print("Step 2: This produces one secondary oocyte with 'XwXw' and one with '0'.")
print("Step 3: After Meiosis II, this leads to the formation of 'XwXw' eggs and '0' eggs.")
print("Result: This single event produces both types of abnormal eggs that explain all of Bridges' exceptional offspring (XXY females and X0 males).")

# Option B: Non-disjunction in Meiosis II
print("\nScenario B: Non-disjunction in Meiosis II")
print("Step 1: Meiosis I is normal. Two secondary oocytes with 'Xw' are formed.")
print("Step 2: In one of these oocytes, the sister chromatids fail to separate in Meiosis II.")
print("Step 3: This lineage produces one 'XwXw' egg and one '0' egg. The other lineage produces normal 'Xw' eggs.")
print("Result: This event can also produce a '0' egg, but non-disjunction in Meiosis I is the more direct event that explains the production of BOTH XX and 0 gametes from a single meiotic division.")

# 6. Conclusion
print("\nConclusion:")
print("The presence of an X0 male requires a '0' egg from the mother. This is caused by non-disjunction of the X chromosomes during female meiosis.")
print("Bridges' observation of both exceptional XXY females and X0 males is best explained by a non-disjunction event during Meiosis I in the female, which creates both the XX and 0 eggs required.")
print("Therefore, the specific event is non-disjunction of the X chromosome in female meiosis I.")
