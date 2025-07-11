def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) for a dry fruit product.
    """
    # 1. Define the given parameters
    net_weight_g = 250.0  # Net weight of the final dry product in grams
    selling_price = 50.0  # Selling price of the product in dollars
    fresh_fruit_moisture_content = 0.20  # 20% moisture in fresh fruits
    
    # Costs associated with the eco-efficiency benchmark for processing a standard unit (assumed to be 1kg)
    benchmark_water_cost = 10.0
    benchmark_energy_cost = 15.0
    
    # Standard unit of mass for the benchmark cost in grams (assumed 1 kg)
    benchmark_mass_g = 1000.0

    # 2. Calculate the required mass of fresh fruit
    # The dry weight is (1 - moisture_content) of the fresh weight.
    # Fresh Weight = Dry Weight / (1 - moisture_content)
    required_fresh_fruit_g = net_weight_g / (1 - fresh_fruit_moisture_content)

    # 3. Calculate the total production cost
    # Total benchmark cost to process 1kg (1000g) of fresh fruit
    total_benchmark_cost_per_kg = benchmark_water_cost + benchmark_energy_cost
    
    # Calculate the proportional cost for our product's required fresh fruit mass
    production_cost = total_benchmark_cost_per_kg * (required_fresh_fruit_g / benchmark_mass_g)

    # 4. Calculate the Sustainable Value Added (SVA)
    # SVA = Revenue - Costs
    sva = selling_price - production_cost
    
    # 5. Print the final equation with all the numbers
    # Using f-string for clear formatting and rounding to 2 decimal places for currency.
    print(f"The Sustainable Value Added (SVA) is calculated as follows:\n")
    print(f"Initial fresh fruit needed for {net_weight_g}g dry product = {required_fresh_fruit_g:.2f}g")
    print(f"Total production cost for {required_fresh_fruit_g:.2f}g of fresh fruit = ${production_cost:.2f}")
    print("\nFinal Equation:")
    print(f"SVA = ${selling_price:.2f} (Revenue) - ${production_cost:.2f} (Total Cost) = ${sva:.2f}")

calculate_sva()