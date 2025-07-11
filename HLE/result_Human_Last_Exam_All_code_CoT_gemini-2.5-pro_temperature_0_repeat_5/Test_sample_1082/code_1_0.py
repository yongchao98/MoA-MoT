def calculate_sva():
    """
    Calculates the Sustainable Value Added based on the problem description.
    """
    # Define the financial values from the problem statement
    selling_price = 50
    water_cost = 10
    energy_cost = 15

    # Calculate the total cost of resources
    total_cost = water_cost + energy_cost

    # Calculate the Sustainable Value Added (SVA)
    sva = selling_price - total_cost

    # Print the explanation and the final equation, showing each number
    print("To find the Sustainable Value Added (SVA), we subtract the total resource costs from the selling price.")
    print("\nThe final equation is:")
    print(f"SVA = Selling Price - (Water Cost + Energy Cost)")
    print(f"SVA = ${selling_price} - (${water_cost} + ${energy_cost})")
    print(f"SVA = ${selling_price} - ${total_cost}")
    print(f"SVA = ${sva}")

# Execute the function to print the result
calculate_sva()