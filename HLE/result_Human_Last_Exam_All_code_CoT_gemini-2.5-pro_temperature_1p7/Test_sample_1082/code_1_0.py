def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided data.
    """
    # 1. Define the variables from the problem description.
    
    # Economic Value Added is the sale price of the final product.
    economic_value = 50  # in dollars

    # Environmental Cost is the sum of the costs for water and energy.
    water_cost = 10      # in dollars
    energy_cost = 15     # in dollars

    # 2. Calculate the total environmental cost.
    total_environmental_cost = water_cost + energy_cost

    # 3. Calculate the Sustainable Value Added (SVA).
    # SVA = Economic Value Added - Environmental Cost
    sva = economic_value - total_environmental_cost

    # 4. Print the result in a clear, step-by-step equation format.
    print("Sustainable Value Added (SVA) Calculation:")
    print("SVA = Economic Value - (Water Cost + Energy Cost)")
    print(f"SVA = ${economic_value} - (${water_cost} + ${energy_cost})")
    print("\nFinal Equation:")
    # The final print statement shows all the numbers in the final equation.
    print(f"{economic_value} - ({water_cost} + {energy_cost}) = {sva}")

# Execute the function
calculate_sustainable_value_added()