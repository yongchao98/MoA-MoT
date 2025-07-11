def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided economic and environmental data.
    """
    # 1. Define the variables from the problem description.
    value_created = 50  # Selling price of the product in $
    water_usage = 100   # Liters of water used
    energy_usage = 40   # kWh of energy used
    water_cost_factor = 10 # Environmental cost factor for water
    energy_cost_factor = 15 # Environmental cost factor for energy

    # 2. Calculate the environmental cost for each resource.
    environmental_cost_water = water_usage * water_cost_factor
    environmental_cost_energy = energy_usage * energy_cost_factor

    # 3. Calculate the total environmental cost.
    total_environmental_cost = environmental_cost_water + environmental_cost_energy

    # 4. Calculate the Sustainable Value Added (SVA).
    sva = value_created - total_environmental_cost

    # 5. Print the breakdown of the calculation as requested.
    print("Sustainable Value Added (SVA) Calculation:")
    print("SVA = Value Created - Total Environmental Cost")
    print("SVA = Value Created - (Water Usage * Water Cost Factor + Energy Usage * Energy Cost Factor)")
    print(f"SVA = {value_created} - ({water_usage} * {water_cost_factor} + {energy_usage} * {energy_cost_factor})")
    print(f"SVA = {value_created} - ({environmental_cost_water} + {environmental_cost_energy})")
    print(f"SVA = {value_created} - {total_environmental_cost}")
    print(f"Final SVA = {sva}")

# Execute the function to see the result.
calculate_sva()