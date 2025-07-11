def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on the problem description.
    """
    # Step 1: Define the known values from the problem statement.
    selling_price = 50  # in dollars
    water_cost = 10     # in dollars for the water used
    energy_cost = 15    # in dollars for the energy used

    # Step 2: Calculate the total environmental cost.
    total_environmental_cost = water_cost + energy_cost

    # Step 3: Calculate the Sustainable Value Added (SVA).
    # SVA is the value created (revenue) minus the environmental costs.
    sustainable_value_added = selling_price - total_environmental_cost

    # Step 4: Print the final equation with all numbers, as requested.
    print("The Sustainable Value Added (SVA) is calculated by subtracting the environmental costs from the product's selling price.")
    print("\nCalculation:")
    print(f"SVA = Selling Price - (Water Cost + Energy Cost)")
    print(f"SVA = ${selling_price} - (${water_cost} + ${energy_cost})")
    print(f"SVA = ${selling_price} - ${total_environmental_cost}")
    print(f"SVA = ${sustainable_value_added}")

if __name__ == '__main__':
    calculate_sva()