def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided financial data.
    """
    # 1. Define the financial values from the problem description.
    # The value created is the selling price of the final product.
    revenue = 50

    # The environmental costs for the process are given directly.
    cost_water = 10
    cost_energy = 15

    # 2. Calculate the total environmental cost.
    total_environmental_cost = cost_water + cost_energy

    # 3. Calculate the Sustainable Value Added (SVA).
    sva = revenue - total_environmental_cost

    # 4. Print the final equation and the result, showing each number.
    print(f"The calculation for Sustainable Value Added (SVA) is:")
    print(f"SVA = Revenue - (Cost of Water + Cost of Energy)")
    print(f"SVA = {revenue} - ({cost_water} + {cost_energy})")
    print(f"SVA = {revenue} - {total_environmental_cost}")
    print(f"The final Sustainable Value Added (SVA) is ${sva}.")

# Execute the function to get the answer.
calculate_sva()
<<<25>>>