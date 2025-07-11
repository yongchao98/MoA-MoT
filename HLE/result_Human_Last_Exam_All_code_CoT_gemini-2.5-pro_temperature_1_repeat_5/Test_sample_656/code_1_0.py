import math

def solve_packing_problem():
    """
    This function solves the energy ball packing problem by iterating through
    possible container configurations to maximize the total stored energy.
    """
    best_config = {
        'energy': 0,
        'n1': 0,
        'n2': 0,
        'L': 0,
        'W': 0,
        'H': 0,
        'SA': 0
    }
    
    # We search for integer multipliers k_l, k_w, k_h for the large ball diameter (4cm)
    # The search range is determined by the surface area constraint. A k value
    # beyond 20 is highly unlikely to produce a valid container.
    search_range = 20
    
    # Assume k_l >= k_w >= k_h to avoid redundant checks of the same box dimensions
    for k_h in range(1, search_range):
        for k_w in range(k_h, search_range):
            for k_l in range(k_w, search_range):
                # Box dimensions based on Simple Cubic packing of 4cm diameter balls
                L = 4.0 * k_l
                W = 4.0 * k_w
                H = 4.0 * k_h
                
                # Surface area of the box
                surface_area = 2 * (L * W + L * H + W * H)
                
                if surface_area > 1050:
                    # Since k_l is increasing, further iterations in this loop will also exceed the SA limit
                    break
                
                # Number of 2-cm radius balls (n2)
                n2 = k_l * k_w * k_h
                
                # Number of 1-cm radius balls (n1)
                # These can fit in the voids at the corners of the 4x4x4 packing cells.
                # The number of such voids is (k_l-1) * (k_w-1) * (k_h-1).
                # We need to handle cases where any k is 1.
                n1 = 0
                if k_l > 1 and k_w > 1 and k_h > 1:
                    n1 = (k_l - 1) * (k_w - 1) * (k_h - 1)
                
                # Total energy
                total_energy = n2 * 20 + n1 * 1
                
                if total_energy > best_config['energy']:
                    best_config = {
                        'energy': total_energy,
                        'n1': n1,
                        'n2': n2,
                        'L': L,
                        'W': W,
                        'H': H,
                        'SA': surface_area
                    }

    # Print the final result
    L = best_config['L']
    W = best_config['W']
    H = best_config['H']
    n1 = best_config['n1']
    n2 = best_config['n2']
    total_energy = best_config['energy']
    
    # Output the description in the format [C]a;b
    container_description = f"box {L}x{W}x{H}"
    final_answer_string = f"[{container_description}]{n1};{n2}"
    print(final_answer_string)
    
    # Output the final equation as requested
    print("\nCalculation of the total potential energy:")
    print(f"Number of 1-cm radius balls (a): {n1}")
    print(f"Number of 2-cm radius balls (b): {n2}")
    print(f"The final equation is: {n2} * 20 MJ + {n1} * 1 MJ = {total_energy} MJ")
    
    # Returning the string for the final answer format
    return final_answer_string


# Execute the function to find and print the solution
final_answer = solve_packing_problem()

# The final answer in the required format
# print(f"\n<<<{final_answer}>>>")