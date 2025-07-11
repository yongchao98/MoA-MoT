def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided metrics.
    """
    # 1. Define all the initial values from the problem statement.
    product_weight_kg = 0.250  # 250g product weight
    selling_price = 50.00       # Selling price of the 250g product

    # Benchmark resource use and cost, assumed for 1kg of output
    benchmark_output_kg = 1.0
    benchmark_water_L = 100.0
    benchmark_energy_kWh = 40.0
    benchmark_water_cost = 10.00
    benchmark_energy_cost = 15.00

    # 2. Scale resources for the 250g (0.25kg) product.
    # Resource use per kg = benchmark use / benchmark output
    # Resource use for product = Resource use per kg * product weight
    water_used_for_product_L = (benchmark_water_L / benchmark_output_kg) * product_weight_kg
    energy_used_for_product_kWh = (benchmark_energy_kWh / benchmark_output_kg) * product_weight_kg

    # 3. Calculate the cost of the resources used for the product.
    # Unit cost = benchmark cost / benchmark resource amount
    # Total cost for product = unit cost * resource used for product
    total_water_cost_for_product = (benchmark_water_cost / benchmark_water_L) * water_used_for_product_L
    total_energy_cost_for_product = (benchmark_energy_cost / benchmark_energy_kWh) * energy_used_for_product_kWh

    # 4. Calculate the total environmental cost.
    total_environmental_cost = total_water_cost_for_product + total_energy_cost_for_product

    # 5. Calculate the Sustainable Value Added (SVA).
    # SVA = Revenue - Monetized Environmental Cost
    sva = selling_price - total_environmental_cost

    # 6. Print the final result in the requested equation format.
    # The format shows: SVA = Selling Price - (Water Cost + Energy Cost)
    print("Sustainable Value Added (SVA) Calculation:")
    print(f"{sva:.2f} = {selling_price:.2f} - ({total_water_cost_for_product:.2f} + {total_energy_cost_for_product:.2f})")

calculate_sva()