def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided economic figures.
    """
    # Step 1: Define the variables from the problem statement.
    selling_price = 50  # The final selling price of the product in dollars.
    water_cost = 10     # The cost of water used in the process in dollars.
    energy_cost = 15    # The cost of energy used in the process in dollars.

    # Step 2: Calculate the total resource cost.
    total_resource_cost = water_cost + energy_cost

    # Step 3: Calculate the Sustainable Value Added (SVA).
    sva = selling_price - total_resource_cost

    # Step 4: Print the final equation and the result.
    print("The Sustainable Value Added (SVA) is the value created minus the environmental resource costs.")
    print("SVA = Selling Price - (Water Cost + Energy Cost)")
    print(f"SVA = ${selling_price} - (${water_cost} + ${energy_cost})")
    print(f"SVA = ${selling_price} - ${total_resource_cost}")
    print(f"The final Sustainable Value Added is: ${sva}")

# Execute the function
calculate_sva()