def calculate_heritability():
    """
    Calculates and compares heritability under two scenarios:
    1. Purely additive genetics (like the rabbits).
    2. Genetics with epistatic interactions (like the other species).
    """
    print("--- Scenario 1: Rabbit Experiment (Purely Additive Genetics) ---")
    # Set variance components. To get H^2 = 0.75, we can set V_A=3 and V_E=1.
    Va_rabbit = 3.0
    Vd_rabbit = 0.0
    Vi_rabbit = 0.0
    Ve_rabbit = 1.0

    # Calculate total genetic and phenotypic variance
    Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
    Vp_rabbit = Vg_rabbit + Ve_rabbit

    # Calculate heritabilities
    H2_rabbit = Vg_rabbit / Vp_rabbit
    h2_rabbit = Va_rabbit / Vp_rabbit

    print(f"Variance Components: Additive(Va)={Va_rabbit}, Dominance(Vd)={Vd_rabbit}, Epistatic(Vi)={Vi_rabbit}, Environmental(Ve)={Ve_rabbit}")
    print(f"Total Genetic Variance (Vg) = {Va_rabbit} + {Vd_rabbit} + {Vi_rabbit} = {Vg_rabbit}")
    print(f"Total Phenotypic Variance (Vp) = {Vg_rabbit} + {Ve_rabbit} = {Vp_rabbit}")
    print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit:.2f}")
    print("Result: With purely additive genetics, H^2 and h^2 are equal.\n")


    print("--- Scenario 2: Other Species (With Epistatic Interactions) ---")
    # Keep other components the same for comparison, but add epistatic variance
    Va_other = 3.0
    Vd_other = 0.0
    Vi_other = 0.8  # Introducing epistatic variance
    Ve_other = 1.0

    # Calculate total genetic and phenotypic variance
    Vg_other = Va_other + Vd_other + Vi_other
    Vp_other = Vg_other + Ve_other

    # Calculate heritabilities
    H2_other = Vg_other / Vp_other
    h2_other = Va_other / Vp_other
    
    print(f"Variance Components: Additive(Va)={Va_other}, Dominance(Vd)={Vd_other}, Epistatic(Vi)={Vi_other}, Environmental(Ve)={Ve_other}")
    print(f"Total Genetic Variance (Vg) = {Va_other} + {Vd_other} + {Vi_other} = {Vg_other}")
    print(f"Total Phenotypic Variance (Vp) = {Vg_other} + {Ve_other} = {Vp_other}")
    print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_other} / {Vp_other} = {H2_other:.2f}")
    print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_other} / {Vp_other} = {h2_other:.2f}")
    print("Result: The presence of epistatic variance (Vi > 0) causes H^2 to be different from (and greater than) h^2.")

# Run the calculation
calculate_heritability()