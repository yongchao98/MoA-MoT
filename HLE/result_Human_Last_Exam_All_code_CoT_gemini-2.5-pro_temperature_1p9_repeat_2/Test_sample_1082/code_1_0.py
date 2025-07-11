def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) based on given revenue and environmental costs.
    """
    # 1. Define the revenue and costs from the problem description.
    revenue = 50  # Revenue from selling the final product in dollars.
    cost_water = 10  # Cost of water used in the process in dollars.
    cost_energy = 15  # Cost of energy used in the process in dollars.

    # 2. Calculate the total environmental cost.
    total_environmental_cost = cost_water + cost_energy

    # 3. Calculate the Sustainable Value Added (SVA).
    sva = revenue - total_environmental_cost

    # 4. Print the final result, showing the equation with the numbers.
    print("To calculate the Sustainable Value Added (SVA), we subtract the total environmental costs from the revenue.")
    print(f"The equation is: Revenue - (Cost of Water + Cost of Energy)")
    print(f"SVA = ${revenue} - (${cost_water} + ${cost_energy}) = ${sva}")


if __name__ == "__main__":
    calculate_sustainable_value_added()