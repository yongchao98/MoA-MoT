def calculate_heritability(V_A, V_D, V_I, V_E):
    """
    Calculates and prints broad-sense (H^2) and narrow-sense (h^2) heritability.
    """
    # Calculate total genetic variance and total phenotypic variance
    V_G = V_A + V_D + V_I
    V_P = V_G + V_E

    # Ensure V_P is not zero to avoid division by zero
    if V_P == 0:
        print("Total phenotypic variance is zero; heritability is undefined.")
        return

    # Calculate heritabilities
    H_squared = V_G / V_P
    h_squared = V_A / V_P

    # Print the results with full equations
    print(f"H^2 = (V_A + V_D + V_I) / (V_A + V_D + V_I + V_E)")
    print(f"H^2 = ({V_A} + {V_D} + {V_I}) / ({V_A} + {V_D} + {V_I} + {V_E}) = {V_G} / {V_P} = {H_squared:.4f}")
    
    print(f"\nh^2 = V_A / (V_A + V_D + V_I + V_E)")
    print(f"h^2 = {V_A} / ({V_A} + {V_D} + {V_I} + {V_E}) = {V_A} / {V_P} = {h_squared:.4f}")
    print("-" * 30)


# --- Scenario 1: The Rabbit Experiment (Entirely Additive) ---
# We set component values to achieve H^2 = 0.75 as given in the problem.
# Let V_A = 3. This means V_G = 3.
# To get H^2 = V_G / (V_G + V_E) = 0.75, then 3 / (3 + V_E) = 0.75 => V_E = 1.
print("Scenario 1: Entirely Additive Model (like the rabbits)")
print("Here, V_D = 0 and V_I = 0, so H^2 and h^2 are equal.\n")
calculate_heritability(V_A=3.0, V_D=0.0, V_I=0.0, V_E=1.0)


# --- Scenario 2: Introducing Epistasis (Testing Choice C) ---
# We keep other values the same but add a non-zero epistatic variance (V_I).
# This represents a situation that would cause a difference between H^2 and h^2.
print("Scenario 2: Introducing Epistatic Interactions")
print("Here, V_I > 0, which causes H^2 and h^2 to differ.\n")
calculate_heritability(V_A=3.0, V_D=0.0, V_I=0.5, V_E=1.0)