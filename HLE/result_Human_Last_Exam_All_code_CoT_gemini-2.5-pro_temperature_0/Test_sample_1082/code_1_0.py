def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) for a product based on batch production data.
    """
    # --- Given variables ---
    product_net_weight_g = 250
    water_use_liters = 100
    energy_use_kwh = 40
    moisture_loss_percentage = 20
    water_cost = 10
    energy_cost = 15
    batch_revenue = 50

    # --- Step 1: Calculate the total environmental cost for the benchmark batch ---
    print("Step 1: Calculate the total environmental cost for the batch.")
    total_environmental_cost = water_cost + energy_cost
    print(f"Total Environmental Cost = Water Cost + Energy Cost = ${water_cost} + ${energy_cost} = ${total_environmental_cost}\n")

    # --- Step 2: Calculate the Sustainable Value Added (SVA) for the entire batch ---
    print("Step 2: Calculate the Sustainable Value Added (SVA) for the batch.")
    sva_batch = batch_revenue - total_environmental_cost
    print(f"SVA for Batch = Batch Revenue - Total Environmental Cost = ${batch_revenue} - ${total_environmental_cost} = ${sva_batch}\n")

    # --- Step 3: Calculate the total weight of the dried fruit batch ---
    # Assumption: 100 liters of water used is the water removed, and 1L of water weighs 1kg.
    water_removed_kg = float(water_use_liters)
    # Assumption: The weight of water removed is 20% of the initial fresh fruit weight.
    initial_fresh_weight_kg = water_removed_kg / (moisture_loss_percentage / 100.0)
    final_batch_weight_kg = initial_fresh_weight_kg - water_removed_kg
    print("Step 3: Calculate the total weight of the final product batch.")
    print(f"Based on the assumption that 100L of water (100 kg) is removed and this represents 20% of the initial fresh weight:")
    print(f"Initial Fresh Weight = {water_removed_kg} kg / {moisture_loss_percentage}% = {initial_fresh_weight_kg} kg")
    print(f"Final Batch Weight = Initial Weight - Water Removed = {initial_fresh_weight_kg} kg - {water_removed_kg} kg = {final_batch_weight_kg} kg\n")

    # --- Step 4: Calculate the SVA for the specific 250g product ---
    product_net_weight_kg = product_net_weight_g / 1000.0
    sva_product = (sva_batch / final_batch_weight_kg) * product_net_weight_kg
    print(f"Step 4: Calculate the SVA for the {product_net_weight_g}g product by scaling the batch SVA.")
    print("Final SVA = (SVA for Batch / Final Batch Weight) * Product Weight")
    # Final equation with all numbers
    print(f"Final SVA = (${sva_batch} / {final_batch_weight_kg} kg) * {product_net_weight_kg} kg\n")

    print("---------------------------------")
    print(f"The final Sustainable Value Added is: ${sva_product:.6f}")
    print("---------------------------------")

calculate_sustainable_value_added()