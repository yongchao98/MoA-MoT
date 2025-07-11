def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided financial data.
    """
    # Step 1: Define the given financial values from the problem description.
    selling_price = 50  # Revenue from selling the product in dollars.
    water_cost = 10     # Cost of water used in dollars.
    energy_cost = 15    # Cost of energy used in dollars.

    # Step 2: Calculate the total environmental cost.
    total_environmental_cost = water_cost + energy_cost

    # Step 3: Calculate the Sustainable Value Added (SVA).
    # SVA = Revenue - Total Environmental Costs
    sustainable_value_added = selling_price - total_environmental_cost

    # Step 4: Print the final equation with each number clearly shown.
    print(f"Calculating Sustainable Value Added (SVA):")
    print(f"SVA = Selling Price - (Cost of Water + Cost of Energy)")
    print(f"SVA = ${selling_price} - (${water_cost} + ${energy_cost})")
    print(f"SVA = ${selling_price} - ${total_environmental_cost}")
    print(f"The final Sustainable Value Added is: ${sustainable_value_added}")

# Execute the function
calculate_sva()