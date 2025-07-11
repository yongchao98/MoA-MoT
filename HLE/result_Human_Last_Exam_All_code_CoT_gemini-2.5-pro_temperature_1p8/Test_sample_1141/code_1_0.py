import math

def solve_cutting_problem():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    into pieces of different shapes.
    """
    
    # Billet dimensions in cm
    billet_length = 16
    billet_width = 11
    billet_height = 4

    # --- Scenario A: Cutting only B2 balls ---
    
    # B2 ball properties
    b2_radius = 2.0
    b2_diameter = b2_radius * 2
    b2_price = 150

    # Since B2 balls take up the full height of the billet, this becomes a 2D packing problem:
    # how many 4cm diameter circles can fit in a 16x11cm rectangle?
    # We can use a simple grid packing for a reliable calculation.
    num_b2_length = math.floor(billet_length / b2_diameter)
    num_b2_width = math.floor(billet_width / b2_diameter)
    num_b2 = num_b2_length * num_b2_width
    value_scenario_a = num_b2 * b2_price

    # --- Scenario B: Cutting only T1 cubes and B1 balls ---

    # T1 cube properties
    t1_side = 1.0
    t1_price = 5

    # B1 ball properties
    b1_diameter = 1.0
    b1_price = 1

    # T1 cubes are more valuable per volume and pack perfectly without wasted space.
    # Therefore, the optimal strategy for this scenario is to fill the billet completely with T1 cubes.
    num_t1_length = math.floor(billet_length / t1_side)
    num_t1_width = math.floor(billet_width / t1_side)
    num_t1_height = math.floor(billet_height / t1_side)
    num_t1 = num_t1_length * num_t1_width * num_t1_height
    value_scenario_b = num_t1 * t1_price
    
    # --- Comparison and Final Result ---

    print("--- Analysis ---")
    print(f"Scenario A (only B2 balls):")
    print(f"Number of B2 balls: {num_b2_length} * {num_b2_width} = {num_b2}")
    print(f"Total Value: {num_b2} * {b2_price} = {value_scenario_a}")
    print("\n")
    print(f"Scenario B (only T1 cubes):")
    print(f"Number of T1 cubes: {num_t1_length} * {num_t1_width} * {num_t1_height} = {num_t1}")
    print(f"Total Value: {num_t1} * {t1_price} = {value_scenario_b}")
    print("\n")
    print("--- Conclusion ---")
    
    if value_scenario_a > value_scenario_b:
        print("Scenario A is optimal.")
        final_value = value_scenario_a
        # This part of the code is for display purposes. It reconstructs the final equation.
        final_equation = f"{num_b2} * {b2_price} = {final_value}"
    else:
        print("Scenario B is optimal.")
        final_value = value_scenario_b
        # This part of the code is for display purposes. It reconstructs the final equation.
        final_equation = f"{num_t1} * {t1_price} = {final_value}"

    print(f"The highest possible value is {final_value}.")
    print("The final calculation is:")
    print(final_equation)


solve_cutting_problem()

# The final answer is the maximum value found.
final_answer = 16 * 11 * 4 * 5
# The 'answer' format below is requested by the prompt.
# <<<3520>>>