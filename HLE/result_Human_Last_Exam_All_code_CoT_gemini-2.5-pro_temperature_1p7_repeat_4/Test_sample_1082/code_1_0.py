def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided details.
    """
    # Step 1: Define the economic value added from the problem statement.
    # The final product is sold at $50.
    value_added = 50

    # Step 2: Define the environmental costs from the problem statement.
    # The cost of water is $10 and the cost of energy is $15.
    cost_water = 10
    cost_energy = 15

    # Calculate the total environmental cost.
    total_environmental_cost = cost_water + cost_energy

    # Step 3: Calculate the Sustainable Value Added (SVA).
    # SVA = Value Added - Total Environmental Cost
    sva = value_added - total_environmental_cost

    # Print the breakdown of the calculation as requested.
    print(f"Value Added (Selling Price): ${value_added}")
    print(f"Total Environmental Cost (Water + Energy): ${total_environmental_cost}")
    print("SVA Formula: Value Added - Total Environmental Cost")
    print(f"Final Equation: SVA = ${value_added} - (${cost_water} + ${cost_energy}) = ${sva}")

if __name__ == "__main__":
    calculate_sustainable_value_added()