# 1. Define the components and heritability formulas

# V_A: Additive genetic variance
# V_D: Dominance genetic variance
# V_I: Epistatic (interaction) genetic variance
# V_E: Environmental variance

def calculate_heritability(V_A, V_D, V_I, V_E):
    """Calculates and prints heritability values based on variance components."""
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E

    # Ensure V_P is not zero to avoid division errors
    if V_P == 0:
        print("Phenotypic variance is zero, heritability is undefined.")
        return

    # Broad-sense heritability (H^2)
    H_squared = V_G / V_P

    # Narrow-sense heritability (h^2)
    h_squared = V_A / V_P

    print(f"Variance Components: V_A={V_A}, V_D={V_D}, V_I={V_I}, V_E={V_E}")
    print(f"Total Genetic Variance (V_G) = V_A + V_D + V_I = {V_A} + {V_D} + {V_I} = {V_G}")
    print(f"Phenotypic Variance (V_P) = V_G + V_E = {V_G} + {V_E} = {V_P}")
    print(f"Broad-Sense Heritability (H^2) = V_G / V_P = {V_G} / {V_P} = {H_squared:.3f}")
    print(f"Narrow-Sense Heritability (h^2) = V_A / V_P = {V_A} / {V_P} = {h_squared:.3f}")
    print("-" * 30)

# 2. Simulate the rabbit scenario from the problem
print("Scenario 1: Rabbit experiment (purely additive genetics)")
# Given H^2 = 0.75. Let's assume V_P = 100 for simplicity.
# V_G = H^2 * V_P = 0.75 * 100 = 75.
# V_E = V_P - V_G = 100 - 75 = 25.
# Since it's purely additive, V_A = V_G = 75, and V_D = 0, V_I = 0.
V_A_rabbit = 75
V_D_rabbit = 0
V_I_rabbit = 0
V_E_rabbit = 25
calculate_heritability(V_A_rabbit, V_D_rabbit, V_I_rabbit, V_E_rabbit)
print("Observation: In a purely additive model, H^2 and h^2 are equal.")

# 3. Test choice C: Introduce epistasis
print("\nScenario 2: Testing Choice C (Presence of epistatic interactions)")
print("Let's add epistatic variance (V_I > 0) to the model.")
# Keep other components the same for a clear comparison
V_A_epistasis = 75
V_D_epistasis = 0
V_I_epistasis = 20 # Introduce epistatic variance
V_E_epistasis = 25
calculate_heritability(V_A_epistasis, V_D_epistasis, V_I_epistasis, V_E_epistasis)
print("Observation: Adding epistatic variance (V_I) causes H^2 and h^2 to differ (H^2 > h^2).")
print("This supports Choice C as a reason for differences in heritability measures.")

# 4. Test choice E: Show impact of dominance variance on h^2
print("\nScenario 3: Testing Choice E's claim (Impact of dominance variance)")
print("Choice E claims dominance variance (V_D) has no impact on h^2. Let's test this.")
# Let's add dominance variance (V_D > 0) to the original rabbit model.
V_A_dominance = 75
V_D_dominance = 15 # Introduce dominance variance
V_I_dominance = 0
V_E_dominance = 25
calculate_heritability(V_A_dominance, V_D_dominance, V_I_dominance, V_E_dominance)
print(f"Observation: The initial h^2 was {0.75:.3f}. With V_D > 0, the new h^2 is different.")
print("This disproves the claim in Choice E that dominance has no impact on h^2; it lowers h^2 by increasing the denominator (V_P).")

print("\nFinal Conclusion:")
print("The primary reason for a difference between broad-sense (H^2) and narrow-sense (h^2) heritability is the presence of non-additive genetic variance (dominance V_D or epistasis V_I).")
print("Choice C correctly identifies that epistatic interactions would cause these measures to differ.")
print("Choice E is incorrect because while dominance variance (V_D) also causes a difference, the choice falsely claims it has no impact on h^2.")