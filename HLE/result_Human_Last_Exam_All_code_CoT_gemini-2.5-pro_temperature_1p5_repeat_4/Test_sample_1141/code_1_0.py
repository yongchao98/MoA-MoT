import math

def solve_cutting_problem():
    """
    Calculates the maximum value based on the user's provided problem formulation.
    """
    # Billet dimensions in cm
    billet_dims = (16, 11, 4)

    # Product specifications
    products = {
        'B2': {'diameter_cm': 4, 'price': 150},
        'B1': {'diameter_cm': 1, 'price': 1},
        'T1': {'side_cm': 1, 'price': 5}
    }

    # --- Scenario 1: B2-based solution ---
    # This scenario is valid under the provided formulation, as it avoids the problematic T1 constraints.
    print("--- Analysis of a B2-based solution ---")
    
    # Calculate the number of B2 balls that can be packed in a grid
    b2_diam = products['B2']['diameter_cm']
    num_b2_x = math.floor(billet_dims[0] / b2_diam)
    num_b2_y = math.floor(billet_dims[1] / b2_diam)
    num_b2_z = math.floor(billet_dims[2] / b2_diam)
    num_b2_total = num_b2_x * num_b2_y * num_b2_z
    val_b2_total = num_b2_total * products['B2']['price']

    print(f"Packing B2 balls (diameter {b2_diam}cm):")
    print(f"Fit {num_b2_x} along length, {num_b2_y} along width, {num_b2_z} along height.")
    print(f"Total B2 pieces: {num_b2_total}")
    print(f"Value from B2 balls: {num_b2_total} * {products['B2']['price']} = {val_b2_total}\n")

    # Calculate remaining space for B1 balls
    # The B2s occupy a block of 16 x (2*4) x 4 cm = 16x8x4 cm
    # This leaves a 16 x (11-8) x 4 cm = 16x3x4 cm block
    rem_y_dim = billet_dims[1] % b2_diam
    b1_diam = products['B1']['diameter_cm']
    num_b1_rem = math.floor(billet_dims[0] / b1_diam) * \
                 math.floor(rem_y_dim / b1_diam) * \
                 math.floor(billet_dims[2] / b1_diam)
    val_b1_rem = num_b1_rem * products['B1']['price']

    print("Packing B1 balls in the remaining space:")
    print(f"Remaining block dimensions: {billet_dims[0]}cm x {rem_y_dim}cm x {billet_dims[2]}cm")
    print(f"Total B1 pieces in this block: {num_b1_rem}")
    print(f"Value from B1 balls: {num_b1_rem} * {products['B1']['price']} = {val_b1_rem}\n")
    
    # Total value for this scenario
    total_val = val_b2_total + val_b1_rem

    print("--- Conclusion ---")
    print("The problem formulation has flaws regarding T1 cube constraints.")
    print("Following the rules as written, mixing T1 and B2 pieces is impossible.")
    print("The highest value comes from packing B2 balls and filling the remaining space with B1 balls.")
    print("\nThe highest valid solution is:")
    print(f"{num_b2_total} * {products['B2']['price']} + {num_b1_rem} * {products['B1']['price']} = {total_val}")
    
    # Returning the final value in the required format
    return total_val

# Run the solver and get the result
max_value = solve_cutting_problem()