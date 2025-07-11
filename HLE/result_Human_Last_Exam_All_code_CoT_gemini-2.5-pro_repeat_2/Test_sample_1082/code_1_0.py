def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided financial data.
    """
    # Step 1: Define the variables from the problem description.
    # Value Added is the revenue from selling the product.
    value_added = 50

    # Environmental costs are the monetized costs of the resources used.
    cost_water = 10
    cost_energy = 15

    # Step 2: Calculate the total environmental cost.
    total_environmental_cost = cost_water + cost_energy

    # Step 3: Calculate the Sustainable Value Added (SVA).
    sva = value_added - total_environmental_cost

    # Step 4: Print the breakdown and the final equation.
    print("Calculating the Sustainable Value Added (SVA):")
    print(f"Value Added (Revenue) = ${value_added}")
    print(f"Total Environmental Cost = Water Cost (${cost_water}) + Energy Cost (${cost_energy}) = ${total_environmental_cost}")
    print("--------------------------------------------------")
    print("SVA = Value Added - Total Environmental Cost")
    # The final equation with each number included
    print(f"SVA = ${value_added} - ${total_environmental_cost}")
    print(f"The final Sustainable Value Added is: ${sva}")

# Execute the function
calculate_sva()