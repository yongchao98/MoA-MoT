def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) for a product based on
    its economic value and the environmental costs of its production.
    """
    # --- 1. Define initial parameters from the problem ---
    product_net_weight_g = 250.0
    fresh_fruit_moisture_content = 0.20
    sale_price = 50.0

    # Benchmark is assumed for processing 1kg of fresh fruit
    benchmark_fresh_fruit_input_kg = 1.0
    # Environmental cost associated with the benchmark resource usage
    benchmark_water_cost = 10.0
    benchmark_energy_cost = 15.0

    # --- 2. Calculate the required weight of fresh fruit ---
    # The final 250g dry product is the solid part (1 - moisture) of the initial fresh fruit.
    fresh_fruit_weight_g = product_net_weight_g / (1 - fresh_fruit_moisture_content)
    # Convert grams to kilograms to match the benchmark unit
    fresh_fruit_weight_kg = fresh_fruit_weight_g / 1000.0

    # --- 3. Calculate the product's specific environmental cost ---
    # The resource consumption is proportional to the weight of fresh fruit processed.
    # We find the ratio of our product's fruit input to the benchmark's fruit input.
    resource_usage_ratio = fresh_fruit_weight_kg / benchmark_fresh_fruit_input_kg

    # Apply this ratio to the benchmark environmental costs
    product_water_cost = resource_usage_ratio * benchmark_water_cost
    product_energy_cost = resource_usage_ratio * benchmark_energy_cost
    total_environmental_cost = product_water_cost + product_energy_cost

    # --- 4. Calculate the Sustainable Value Added (SVA) ---
    sva = sale_price - total_environmental_cost

    # --- 5. Output the results including the final equation ---
    print("Calculation of Sustainable Value Added (SVA):")
    print("-" * 40)
    print(f"The economic value (sale price) of the product is ${sale_price:.2f}.")
    print(f"To produce {product_net_weight_g:.1f}g of dry fruit, {fresh_fruit_weight_g:.1f}g of fresh fruit is needed.")

    print("\nThe environmental cost is calculated proportionally to the benchmark:")
    print(f"Environmental Cost of Water = ({fresh_fruit_weight_kg:.4f} kg / {benchmark_fresh_fruit_input_kg:.1f} kg) * ${benchmark_water_cost:.2f} = ${product_water_cost:.4f}")
    print(f"Environmental Cost of Energy = ({fresh_fruit_weight_kg:.4f} kg / {benchmark_fresh_fruit_input_kg:.1f} kg) * ${benchmark_energy_cost:.2f} = ${product_energy_cost:.4f}")
    
    print("\nThe final SVA equation is:")
    print(f"SVA = Sale Price - (Water Cost + Energy Cost)")
    # The instruction is to output each number in the final equation.
    print(f"SVA = {sale_price} - ({product_water_cost} + {product_energy_cost})")
    print(f"SVA = {sale_price} - {total_environmental_cost}")
    print(f"SVA = {sva}")


# Run the calculation and print the results
calculate_sustainable_value_added()