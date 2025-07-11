def calculate_heritability():
    """
    This script explains and demonstrates the concepts of heritability to solve the problem.
    """
    print("--- Understanding Heritability Formulas ---")
    print("Broad-Sense Heritability (H^2) = Vg / Vp = (Va + Vd + Vi) / (Va + Vd + Vi + Ve)")
    print("Narrow-Sense Heritability (h^2) = Va / Vp = Va / (Va + Vd + Vi + Ve)")
    print("Where:")
    print("  Va = Additive Genetic Variance")
    print("  Vd = Dominance Genetic Variance")
    print("  Vi = Epistatic (Interaction) Variance")
    print("  Ve = Environmental Variance")
    print("  Vg = Total Genetic Variance (Va + Vd + Vi)")
    print("  Vp = Total Phenotypic Variance (Vg + Ve)")
    print("-" * 40)

    # Step 1: Model the Rabbit Experiment
    print("--- Rabbit Scenario Analysis ---")
    print("Given: H^2 = 0.75 and genetic variance is entirely additive.")
    print("This means Vd = 0 and Vi = 0, so H^2 should equal h^2.")
    # We can set values where Va = 3 * Ve to get H^2 = 0.75
    Va_rabbit = 3.0
    Vd_rabbit = 0.0
    Vi_rabbit = 0.0
    Ve_rabbit = 1.0
    Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
    Vp_rabbit = Vg_rabbit + Ve_rabbit
    H2_rabbit = Vg_rabbit / Vp_rabbit
    h2_rabbit = Va_rabbit / Vp_rabbit
    print(f"Rabbit Va={Va_rabbit}, Vd={Vd_rabbit}, Vi={Vi_rabbit}, Ve={Ve_rabbit}")
    print(f"Rabbit Vg={Vg_rabbit}, Vp={Vp_rabbit}")
    print(f"Rabbit H^2 = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit:.2f}")
    print(f"Rabbit h^2 = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit:.2f}")
    print("As expected, H^2 and h^2 are equal.\n")
    print("-" * 40)

    # Step 2: Analyze what causes H^2 and h^2 to differ
    print("--- Evaluating Choice C: Presence of Epistatic Interactions ---")
    print("Let's model another species where epistatic variance (Vi) is present.")
    print("We will add Vi and see how it affects the heritability measures.")
    
    # We use the rabbit's values as a base and add epistatic variance
    Va_species = 3.0
    Vd_species = 0.0
    Vi_species = 1.0  # Introducing epistatic variance
    Ve_species = 1.0

    Vg_species = Va_species + Vd_species + Vi_species
    Vp_species = Vg_species + Ve_species
    
    H2_species = Vg_species / Vp_species
    h2_species = Va_species / Vp_species
    
    print(f"New Species Scenario (with Vi > 0):")
    print(f"  Va = {Va_species}, Vd = {Vd_species}, Vi = {Vi_species}, Ve = {Ve_species}")
    print("\nCalculating Total Genetic Variance (Vg):")
    print(f"  Vg = Va + Vd + Vi = {Va_species} + {Vd_species} + {Vi_species} = {Vg_species}")
    
    print("\nCalculating Total Phenotypic Variance (Vp):")
    print(f"  Vp = Vg + Ve = {Vg_species} + {Ve_species} = {Vp_species}")

    print("\nCalculating Heritability:")
    print(f"  Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_species} / {Vp_species} = {H2_species:.2f}")
    print(f"  Narrow-Sense Heritability (h^2) = Va / Vp = {Va_species} / {Vp_species} = {h2_species:.2f}")

    print("\nConclusion:")
    print(f"The presence of epistatic variance (Vi = {Vi_species}) caused H^2 ({H2_species:.2f}) and h^2 ({h2_species:.2f}) to be different.")
    print("This confirms that Choice C is a valid explanation for differences in heritability measures.")
    print("-" * 40)
    print("Note on Choice E: Introducing dominance variance (Vd) would also make H^2 and h^2 differ, but the statement claims it has 'no impact on h^2', which is false as Vd increases the denominator (Vp), thus lowering h^2.")

# Execute the function to print the analysis
calculate_heritability()