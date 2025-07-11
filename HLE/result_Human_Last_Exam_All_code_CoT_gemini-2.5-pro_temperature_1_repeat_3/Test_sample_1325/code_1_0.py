def calculate_heritability():
    """
    Calculates and explains broad-sense (H^2) and narrow-sense (h^2) heritability
    under different genetic variance scenarios.
    """

    print("--- Scenario 1: Rabbits (Purely Additive Genetics) ---")
    # Define variance components for the rabbit population
    # Let's assume a total phenotypic variance (Vp) of 100 for simplicity.
    # Given H^2 = 0.75, Vg must be 75.
    # Since genetics are 'entirely additive', Vd and Vi are 0.
    Vp_rabbit = 100
    Va_rabbit = 75  # Additive variance
    Vd_rabbit = 0   # Dominance variance
    Vi_rabbit = 0   # Epistatic variance
    Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
    Ve_rabbit = Vp_rabbit - Vg_rabbit

    # Calculate heritabilities for rabbits
    H2_rabbit = Vg_rabbit / Vp_rabbit
    h2_rabbit = Va_rabbit / Vp_rabbit

    print(f"The components of variance are:")
    print(f"Additive (Va) = {Va_rabbit}, Dominance (Vd) = {Vd_rabbit}, Epistatic (Vi) = {Vi_rabbit}")
    print(f"Total Genetic (Vg = Va + Vd + Vi) = {Vg_rabbit}")
    print(f"Environmental (Ve) = {Ve_rabbit}")
    print(f"Total Phenotypic (Vp = Vg + Ve) = {Vp_rabbit}\n")

    print("Calculating heritability for rabbits:")
    print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit:.2f}")
    print("As you can see, when genetic variance is purely additive, H^2 and h^2 are equal.\n")

    print("--- Scenario 2: Other Species (with Epistasis) ---")
    # Let's introduce epistatic variance (Vi) while keeping total variance the same
    # This means some of the genetic variance is now non-additive
    Vp_other = 100
    Vi_other = 15   # Introduce epistatic variance
    Vd_other = 0    # Assume no dominance variance for this example
    Va_other = 60   # Additive variance is now smaller
    Vg_other = Va_other + Vd_other + Vi_other
    Ve_other = Vp_other - Vg_other

    # Calculate heritabilities for the other species
    H2_other = Vg_other / Vp_other
    h2_other = Va_other / Vp_other

    print(f"The components of variance are:")
    print(f"Additive (Va) = {Va_other}, Dominance (Vd) = {Vd_other}, Epistatic (Vi) = {Vi_other}")
    print(f"Total Genetic (Vg = Va + Vd + Vi) = {Vg_other}")
    print(f"Environmental (Ve) = {Ve_other}")
    print(f"Total Phenotypic (Vp = Vg + Ve) = {Vp_other}\n")

    print("Calculating heritability for the other species:")
    print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_other} / {Vp_other} = {H2_other:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_other} / {Vp_other} = {h2_other:.2f}")
    print("\nConclusion:")
    print("The presence of non-additive genetic variance (like epistasis, Vi) causes H^2 and h^2 to differ.")
    print("This demonstrates that the presence of epistatic interactions is a valid reason for the heritability measures to be different, supporting answer choice C.")

calculate_heritability()
<<<C>>>