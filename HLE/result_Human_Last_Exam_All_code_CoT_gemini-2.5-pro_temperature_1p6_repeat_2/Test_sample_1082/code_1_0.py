def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on provided eco-efficiency data.
    """
    # Step 1: Define the variables from the problem description
    revenue = 50.0  # Selling price of the final product in dollars
    
    # Benchmark (optimal) resource usage
    benchmark_water_use = 100.0  # liters
    benchmark_energy_use = 40.0  # kWh
    
    # Actual resource usage
    actual_water_use = 10.0  # liters
    actual_energy_use = 15.0  # kWh

    # Step 2: Calculate the benchmark eco-efficiency for each resource
    # This is the value generated per unit of resource in the optimal scenario
    benchmark_eco_efficiency_water = revenue / benchmark_water_use
    benchmark_eco_efficiency_energy = revenue / benchmark_energy_use

    # Step 3: Calculate the environmental cost of the actual resources used
    # This represents the value that should have been created by the resources, per the benchmark
    environmental_cost_water = actual_water_use * benchmark_eco_efficiency_water
    environmental_cost_energy = actual_energy_use * benchmark_eco_efficiency_energy
    total_environmental_cost = environmental_cost_water + environmental_cost_energy

    # Step 4: Calculate the Sustainable Value Added (SVA)
    sva = revenue - total_environmental_cost

    # Step 5: Print the breakdown and the final equation
    print("--- Sustainable Value Added (SVA) Calculation ---")
    print(f"Revenue (Value Added): ${revenue:.2f}")
    print(f"Total Environmental Cost: ${total_environmental_cost:.2f}")
    print("\nFinal Equation:")
    print(f"SVA = Revenue - (Actual Water Use * (Revenue / Benchmark Water Use) + Actual Energy Use * (Revenue / Benchmark Energy Use))")
    print(f"SVA = {revenue} - (({actual_water_use} * ({revenue} / {benchmark_water_use})) + ({actual_energy_use} * ({revenue} / {benchmark_energy_use})))")
    print(f"SVA = {revenue} - ({environmental_cost_water:.2f} + {environmental_cost_energy:.2f})")
    print(f"SVA = {revenue} - {total_environmental_cost:.2f}")
    print(f"Calculated SVA: ${sva:.2f}")

# Execute the function
if __name__ == "__main__":
    calculate_sva()
    print("\n<<<26.25>>>")