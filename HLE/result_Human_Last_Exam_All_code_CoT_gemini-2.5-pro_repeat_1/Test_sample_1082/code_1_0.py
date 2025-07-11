def calculate_sva():
    """
    Calculates the Sustainable Value Added based on the provided project details.
    """
    # Step 1: Define the variables from the problem statement.
    value_added = 50  # Selling price of the final product in dollars
    water_usage = 100  # Liters of water used
    energy_usage = 40  # kWh of energy used
    water_cost_per_unit = 10  # Cost in dollars per liter of water
    energy_cost_per_unit = 15  # Cost in dollars per kWh of energy

    # Step 2: Calculate the total environmental impact cost.
    total_water_cost = water_usage * water_cost_per_unit
    total_energy_cost = energy_usage * energy_cost_per_unit
    environmental_impact_cost = total_water_cost + total_energy_cost

    # Step 3: Calculate the Sustainable Value Added (SVA).
    sva = value_added - environmental_impact_cost

    # Step 4: Print the final equation with all numbers and the result.
    print("Sustainable Value Added (SVA) Calculation:")
    print(f"SVA = Value Added - (Cost of Water + Cost of Energy)")
    print(f"SVA = {value_added} - (({water_usage} * {water_cost_per_unit}) + ({energy_usage} * {energy_cost_per_unit}))")
    print(f"SVA = {value_added} - ({total_water_cost} + {total_energy_cost})")
    print(f"SVA = {value_added} - {environmental_impact_cost}")
    print(f"The final Sustainable Value Added is: ${sva}")

calculate_sva()
<<<-1550>>>