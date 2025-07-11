def calculate_heritability():
    """
    Calculates and explains broad-sense (H^2) and narrow-sense (h^2) heritability
    under different genetic variance scenarios to solve the problem.
    """

    print("--- Understanding Heritability: H^2 vs. h^2 ---\n")
    print("Phenotypic Variance (Vp) = Genetic Variance (Vg) + Environmental Variance (Ve)")
    print("Genetic Variance (Vg) = Additive (Va) + Dominance (Vd) + Epistatic (Vi)\n")
    print("Broad-sense Heritability (H^2) = Vg / Vp")
    print("Narrow-sense Heritability (h^2) = Va / Vp\n")
    print("The key difference is that H^2 includes all genetic variance (Va, Vd, Vi),")
    print("while h^2 only includes the additive portion (Va).\n")

    # --- Scenario 1: Rabbits (Purely Additive Genetics) ---
    print("--- Scenario 1: Rabbit Population (Entirely Additive) ---")
    print("The problem states H^2 = 0.75 and the genetic variance is entirely additive.")
    print("This means Dominance (Vd) and Epistatic (Vi) variances are 0.\n")

    # Let's set component values that result in H^2 = 0.75
    Va_rabbit = 75
    Vd_rabbit = 0
    Vi_rabbit = 0
    Ve_rabbit = 25

    Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
    Vp_rabbit = Vg_rabbit + Ve_rabbit
    H2_rabbit = Vg_rabbit / Vp_rabbit
    h2_rabbit = Va_rabbit / Vp_rabbit

    print(f"Let's assume the following variance components for the rabbits:")
    print(f"Additive Variance (Va) = {Va_rabbit}")
    print(f"Dominance Variance (Vd) = {Vd_rabbit}")
    print(f"Epistatic Variance (Vi) = {Vi_rabbit}")
    print(f"Environmental Variance (Ve) = {Ve_rabbit}\n")

    print("Calculations for the rabbit population:")
    print(f"Total Genetic Variance (Vg) = {Va_rabbit} + {Vd_rabbit} + {Vi_rabbit} = {Vg_rabbit}")
    print(f"Total Phenotypic Variance (Vp) = {Vg_rabbit} + {Ve_rabbit} = {Vp_rabbit}")
    print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit:.2f}\n")
    print("As you can see, when genetic variance is purely additive, H^2 and h^2 are equal.\n")


    # --- Scenario 2: Other Species (with Epistasis) ---
    print("--- Scenario 2: Another Species with Epistatic Interactions (Choice C) ---")
    print("Now, let's analyze what happens if a species has non-additive variance, like epistasis (Vi > 0).\n")

    # We introduce epistatic variance
    Va_other = 75
    Vd_other = 0  # We can keep this 0 to isolate the effect of epistasis
    Vi_other = 15 # Introducing epistatic variance
    Ve_other = 25

    Vg_other = Va_other + Vd_other + Vi_other
    Vp_other = Vg_other + Ve_other
    H2_other = Vg_other / Vp_other
    h2_other = Va_other / Vp_other

    print(f"Let's assume the following variance components for this other species:")
    print(f"Additive Variance (Va) = {Va_other}")
    print(f"Dominance Variance (Vd) = {Vd_other}")
    print(f"Epistatic Variance (Vi) = {Vi_other}  <-- This component is now non-zero")
    print(f"Environmental Variance (Ve) = {Ve_other}\n")

    print("Calculations for the other species:")
    print(f"Total Genetic Variance (Vg) = {Va_other} + {Vd_other} + {Vi_other} = {Vg_other}")
    print(f"Total Phenotypic Variance (Vp) = {Vg_other} + {Ve_other} = {Vp_other}")
    print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_other} / {Vp_other} = {H2_other:.3f}")
    print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_other} / {Vp_other} = {h2_other:.3f}\n")

    print("Conclusion:")
    print(f"In this scenario, H^2 ({H2_other:.3f}) is greater than h^2 ({h2_other:.3f}).")
    print("This difference is caused by the presence of the epistatic variance (Vi).")
    print("Therefore, the presence of epistatic interactions (or dominance variance) causes the two heritability measures to differ.")
    print("This confirms that Choice C is the correct answer.")

# Execute the function
calculate_heritability()