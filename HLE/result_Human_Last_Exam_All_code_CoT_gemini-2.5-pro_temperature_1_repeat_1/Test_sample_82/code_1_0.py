# Plan:
# 1. Define total phenotypic variance (Vp) for a hypothetical population.
# 2. Use the given broad-sense heritability (H2 = 0.5) to calculate total genetic (Vg) and environmental (Ve) variance.
# 3. Assume a plausible breakdown of the total genetic variance (Vg) into its additive (Va), dominance (Vd), and epistatic (Vi) components.
# 4. Calculate the narrow-sense heritability (h2), which is the theoretical maximum variance a standard PGS can explain.
# 5. Print the breakdown and results to illustrate why statements A and C are true.

# Let's assume a total phenotypic variance for a trait
Vp = 100

# Given: Broad-sense heritability (H2) is 0.5
H2 = 0.5

# Step 1: Calculate the total genetic and environmental variance components
Vg = H2 * Vp
Ve = Vp - Vg

# Step 2: Assume a plausible breakdown of the total genetic variance (Vg).
# For most complex traits, not all genetic variance is purely additive.
# We will assign some variance to dominance and epistasis.
# Vg = Va + Vd + Vi
Va = 35  # Additive variance
Vd = 10  # Dominance variance
Vi = 5   # Epistatic (gene-gene interaction) variance

# Sanity check: Ensure the components sum to the total genetic variance
if (Va + Vd + Vi) != Vg:
    print(f"Error: Assumed genetic components do not sum to Vg. {Va}+{Vd}+{Vi} != {Vg}")
else:
    # Step 3: Calculate the narrow-sense heritability (h2)
    # This represents the proportion of variance explained by additive effects,
    # which is the maximum explanatory power of a standard PGS.
    h2 = Va / Vp

    # Step 4: Print the full breakdown of the variance equation
    print("Phenotypic Variance Equation: Vp = Va + Vd + Vi + Ve")
    print(f"Values: {Vp} = {Va} + {Vd} + {Vi} + {Ve}")
    print("-" * 30)

    # Print the calculated heritabilities
    print(f"Broad-Sense Heritability (H2) = Vg / Vp = ({Va} + {Vd} + {Vi}) / {Vp} = {Vg / Vp}")
    print(f"Narrow-Sense Heritability (h2) = Va / Vp = {Va} / {Vp} = {h2}")
    print("-" * 30)

    # Explain the implications for the Polygenic Score (PGS)
    print("Conclusion based on this model:")
    print(f"1. A perfect PGS would explain h2 = {h2:.2f} (or {h2*100}%) of the variance.")
    print(f"2. This value ({h2:.2f}) is not more than H2 (0.5). (Statement A is True)")
    print(f"3. The PGS does not approach an explained variance of 0.5 because dominance and epistasis account for part of the genetic variance. (Statement C is True)")
