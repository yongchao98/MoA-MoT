import math

def calculate_cost_of_gene_flow():
    """
    This function demonstrates the calculation of the cost of gene flow in yeast
    based on the principles outlined in the correct answer choice.

    The cost is quantified by the selection coefficient (s), which measures the
    fitness reduction of hybrids compared to their parental lines.
    """

    # --- Hypothetical Fitness Data ---
    # Fitness measurements for a parental line (no gene flow)
    parent_growth_rate = 0.55  # divisions per hour
    parent_biomass = 12.0    # optical density at 600nm (OD600)
    parent_mating_efficiency = 0.95 # 95% success rate

    # Fitness measurements for a hybrid line (result of gene flow)
    hybrid_growth_rate = 0.51
    hybrid_biomass = 11.2
    hybrid_mating_efficiency = 0.78

    print("To measure the cost of gene flow, we calculate the selection coefficient (s) for hybrids relative to their parents.")
    print("The formula is: s = 1 - (Hybrid_Fitness / Parent_Fitness)\n")

    # --- 1. Calculation for Growth Rate ---
    w_growth = hybrid_growth_rate / parent_growth_rate
    s_growth = 1 - w_growth
    print("Cost based on Growth Rate:")
    print(f"  Parent Fitness = {parent_growth_rate}")
    print(f"  Hybrid Fitness = {hybrid_growth_rate}")
    print(f"  Equation: s = 1 - ({hybrid_growth_rate} / {parent_growth_rate})")
    print(f"  Result: s = {s_growth:.4f} (a {s_growth:.2%} fitness cost)\n")

    # --- 2. Calculation for Biomass Production ---
    w_biomass = hybrid_biomass / parent_biomass
    s_biomass = 1 - w_biomass
    print("Cost based on Biomass Production:")
    print(f"  Parent Fitness = {parent_biomass}")
    print(f"  Hybrid Fitness = {hybrid_biomass}")
    print(f"  Equation: s = 1 - ({hybrid_biomass} / {parent_biomass})")
    print(f"  Result: s = {s_biomass:.4f} (a {s_biomass:.2%} fitness cost)\n")

    # --- 3. Calculation for Mating Efficiency ---
    w_mating = hybrid_mating_efficiency / parent_mating_efficiency
    s_mating = 1 - w_mating
    print("Cost based on Mating Efficiency:")
    print(f"  Parent Fitness = {parent_mating_efficiency}")
    print(f"  Hybrid Fitness = {hybrid_mating_efficiency}")
    print(f"  Equation: s = 1 - ({hybrid_mating_efficiency} / {parent_mating_efficiency})")
    print(f"  Result: s = {s_mating:.4f} (a {s_mating:.2%} fitness cost)\n")

calculate_cost_of_gene_flow()
<<<E>>>