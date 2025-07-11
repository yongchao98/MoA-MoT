def calculate_sustainable_value_added():
    """
    Calculates the sustainable value added based on provided financial data.
    """
    # Step 1: Define the revenue from selling the product.
    selling_price = 50

    # Step 2: Define the environmental costs.
    # The cost of 100 liters of water is $10.
    water_cost = 10
    # The cost of 40 kWh of energy is $15.
    energy_cost = 15

    # Step 3: Calculate the total environmental cost.
    total_environmental_cost = water_cost + energy_cost

    # Step 4: Calculate the Sustainable Value Added (SVA).
    sustainable_value_added = selling_price - total_environmental_cost

    # Step 5: Print the calculation process clearly.
    # The problem requests to output each number in the final equation.
    print(f"Sustainable Value Added = Selling Price - (Water Cost + Energy Cost)")
    print(f"Calculation: {selling_price} - ({water_cost} + {energy_cost})")
    print(f"Final Equation: {selling_price} - {total_environmental_cost} = {sustainable_value_added}")
    print(f"The Sustainable Value Added is: ${sustainable_value_added}")

calculate_sustainable_value_added()