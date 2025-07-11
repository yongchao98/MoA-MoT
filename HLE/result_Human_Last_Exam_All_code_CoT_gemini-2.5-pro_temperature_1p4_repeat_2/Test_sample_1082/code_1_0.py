import sys

def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on the problem description.

    The calculation follows these steps:
    1.  Determines the amount of fresh fruit needed to produce the final dry product.
    2.  Calculates the cost of water and energy to process that amount of fresh fruit,
        based on the provided benchmark costs.
    3.  Calculates the SVA by subtracting the total resource costs from the sale price.
    """

    # --- Parameters from the problem ---
    net_weight_dry_g = 250.0
    sale_price = 50.0
    fresh_fruit_moisture = 0.20

    # Benchmark assumes costs are for processing 1kg of fresh fruit
    benchmark_unit_kg = 1.0
    benchmark_water_cost = 10.0
    benchmark_energy_cost = 15.0

    # --- Step 1: Calculate the required weight of fresh fruit ---
    # The dry matter is 100% - 20% moisture = 80% of the fresh fruit's weight.
    dry_matter_ratio = 1.0 - fresh_fruit_moisture
    required_fresh_fruit_g = net_weight_dry_g / dry_matter_ratio
    required_fresh_fruit_kg = required_fresh_fruit_g / 1000.0

    # --- Step 2: Calculate the total cost of resources ---
    # Cost to process our required amount of fresh fruit is proportional to the benchmark.
    total_water_cost = (required_fresh_fruit_kg / benchmark_unit_kg) * benchmark_water_cost
    total_energy_cost = (required_fresh_fruit_kg / benchmark_unit_kg) * benchmark_energy_cost

    # --- Step 3: Calculate the final SVA ---
    sustainable_value_added = sale_price - (total_water_cost + total_energy_cost)

    # --- Step 4: Output the final equation and the result ---
    print("The final equation for Sustainable Value Added (SVA) is:")
    # Using format specifiers to ensure clarity in the output equation
    # SVA = Revenue - (Water Cost + Energy Cost)
    print(f"{sustainable_value_added:.4f} = {sale_price:.2f} - ({total_water_cost:.4f} + {total_energy_cost:.4f})")
    
    # Final answer in the requested format
    sys.stdout.write(f"\n<<<{sustainable_value_added:.4f}>>>\n")

if __name__ == '__main__':
    calculate_sva()