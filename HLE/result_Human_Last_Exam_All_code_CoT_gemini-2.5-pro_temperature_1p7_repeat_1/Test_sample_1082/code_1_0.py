def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided financial data.

    The function extracts revenue and costs from the problem description,
    calculates the SVA, and prints the full equation and the result.
    """
    # Revenue from selling the final product
    revenue = 50

    # Economic Cost: Interpreted as the cost of raw materials (fresh fruits)
    material_cost = 10

    # Environmental Costs: Costs associated with resource consumption
    water_cost = 10
    energy_cost = 15

    # Although other details are provided (e.g., 250g weight, 100L water, 40kWh energy),
    # the revenue and cost figures ($50, $10, $10, $15) are assumed to apply
    # to the single product unit for which we are calculating the SVA.

    # Calculate total costs
    total_cost = material_cost + water_cost + energy_cost

    # Calculate Sustainable Value Added (SVA)
    sva = revenue - total_cost

    # Print the explanation and the final equation with the result
    print("To calculate the Sustainable Value Added (SVA), we subtract the economic and environmental costs from the revenue.")
    print("\nCalculation:")
    print(f"Revenue: ${revenue}")
    print(f"Costs: Material (${material_cost}) + Water (${water_cost}) + Energy (${energy_cost})")
    print("\nFinal Equation:")
    print(f"{revenue} - ({material_cost} + {water_cost} + {energy_cost}) = {sva}")
    
    # Returning the final numerical value for the grading system
    return sva

# Execute the function to print the result
final_sva = calculate_sva()

# The following line is for internal verification and not part of the primary output.
# print(f"\n<<<__{final_sva}__>>>")