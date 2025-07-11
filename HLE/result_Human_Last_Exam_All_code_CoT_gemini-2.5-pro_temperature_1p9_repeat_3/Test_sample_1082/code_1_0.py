def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided details.
    
    The SVA is calculated as the revenue generated from a product minus the
    environmental costs incurred during its production.
    """
    
    # Define the given financial values from the problem description
    revenue = 50
    water_cost = 10
    energy_cost = 15

    # Calculate the total environmental cost by summing the individual resource costs
    total_environmental_cost = water_cost + energy_cost

    # Calculate the Sustainable Value Added (SVA)
    sustainable_value_added = revenue - total_environmental_cost

    # --- Output the results ---
    print("Step 1: Identify Revenue and Environmental Costs")
    print(f"Revenue from the final product: ${revenue}")
    print(f"Cost of water resource: ${water_cost}")
    print(f"Cost of energy resource: ${energy_cost}")
    print("-" * 40)

    print("Step 2: Calculate the Final Equation for Sustainable Value Added (SVA)")
    print("Formula: SVA = Revenue - (Water Cost + Energy Cost)")
    print(f"Final Equation: SVA = ${revenue} - (${water_cost} + ${energy_cost})")
    print(f"SVA = ${revenue} - ${total_environmental_cost}")
    print(f"The Sustainable Value Added is: ${sustainable_value_added}")

# Execute the calculation
calculate_sva()