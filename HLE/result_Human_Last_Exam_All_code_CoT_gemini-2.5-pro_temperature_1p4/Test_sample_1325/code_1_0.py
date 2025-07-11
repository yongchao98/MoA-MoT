def calculate_and_print_heritability(name, V_A, V_D, V_I, V_E):
    """Calculates and prints heritability components and values."""
    V_G = V_A + V_D + V_I
    V_P = V_A + V_D + V_I + V_E

    if V_P == 0:
        print(f"\nCannot calculate heritability for {name}: Total Phenotypic Variance is zero.")
        return

    H2 = V_G / V_P
    h2 = V_A / V_P

    print(f"\n--- Analysis for {name} ---")
    print(f"Variance Components: V_A={V_A}, V_D={V_D}, V_I={V_I}, V_E={V_E}")
    print(f"Total Genetic Variance (V_G) = {V_A} + {V_D} + {V_I} = {V_G}")
    print(f"Total Phenotypic Variance (V_P) = {V_A} + {V_D} + {V_I} + {V_E} = {V_P}")
    print("\nBroad-Sense Heritability (H^2):")
    print(f"H^2 = V_G / V_P = {V_G} / {V_P} = {H2:.2f}")

    print("\nNarrow-Sense Heritability (h^2):")
    print(f"h^2 = V_A / V_P = {V_A} / {V_P} = {h2:.2f}")

    if H2 == h2:
        print("\nObservation: H^2 and h^2 are equal because there is no non-additive genetic variance (V_D and V_I are zero).")
    else:
        print("\nObservation: H^2 and h^2 differ. This difference is caused by non-additive genetic variance (V_D or V_I).")
        print(f"Note how Dominance Variance (V_D = {V_D}) contributes to the numerator of H^2 but not h^2, causing the divergence.")

# Case 1: Rabbits with entirely additive genetic variance
# We set V_A=3 and V_E=1 to get the H^2 = 0.75 described in the problem.
# V_D and V_I are zero by definition.
calculate_and_print_heritability(name="Rabbits (Purely Additive Model)", V_A=3, V_D=0, V_I=0, V_E=1)

# Case 2: Another species with dominance variance (Illustrating Option E)
# We add dominance variance (V_D=1) to show its effect.
calculate_and_print_heritability(name="Another Species (with Dominance)", V_A=3, V_D=1, V_I=0, V_E=1)
