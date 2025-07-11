def calculate_sva():
    """
    Calculates the Sustainable Value Added based on the provided financial figures.
    """
    # Step 1 & 2: Define the key financial figures from the problem description.
    selling_price = 50
    cost_water = 10
    cost_energy = 15

    # Step 3: Calculate the total cost and the Sustainable Value Added (SVA).
    total_cost = cost_water + cost_energy
    sva = selling_price - total_cost

    # Step 5: Print the final equation with each number and the result.
    print("The Sustainable Value Added (SVA) is the value generated minus the costs of the resources used.")
    print("The final equation is:")
    print(f"{selling_price} - ({cost_water} + {cost_energy}) = {sva}")
    print(f"\nThe Sustainable Value Added for the product is ${sva}.")

# Execute the function to display the result.
calculate_sva()