import math

def solve_pioneer_probe_packing():
    """
    Solves the energy ball packing problem for the Pioneer probe.

    This script finds the optimal container design and ball configuration to maximize
    the total stored energy, subject to a surface area constraint of 1,050 cm^2.

    The strategy is based on a simplified cubic lattice packing model, which prioritizes
    the more energy-dense 2-cm radius balls.
    """
    max_surface_area = 1050
    
    best_n2 = 0
    best_dimensions = (0, 0, 0)
    best_abc = (0, 0, 0)

    # We assume container dimensions L, W, H are multiples of 4 (diameter of 2-cm ball)
    # to achieve the most efficient lattice packing with no wasted edge space.
    # Let L=4a, W=4b, H=4c. The number of 2-cm balls is n2 = a*b*c.
    # The surface area is SA = 32 * (ab + bc + ca).
    # We search for integers a, b, c that maximize n2 subject to SA <= 1050.
    # This implies ab + bc + ca <= 1050 / 32 = 32.8125.
    
    # A search limit of 10 for a, b, c is more than sufficient.
    # E.g., if a=b=c, 3*a^2 < 33 => a^2 < 11 => a <= 3.
    for a in range(1, 10):
        for b in range(a, 10):  # Start from 'a' to avoid redundant permutations
            for c in range(b, 10): # Start from 'b' for the same reason
                
                sa_factor = a * b + b * c + c * a
                current_sa = 32 * sa_factor
                
                if current_sa <= max_surface_area:
                    num_2cm_balls = a * b * c
                    if num_2cm_balls > best_n2:
                        best_n2 = num_2cm_balls
                        best_dimensions = (4 * a, 4 * b, 4 * c)
                        best_abc = (a, b, c)

    # With the optimal dimensions found, the container is perfectly filled with 2-cm balls
    # in our model, leaving no room for the less energy-dense 1-cm balls.
    num_1cm_balls = 0
    num_2cm_balls = best_n2
    l, w, h = best_dimensions
    
    # Construct the final descriptive string and calculate final values.
    container_description = f"box {l}x{w}x{h}"
    final_surface_area = 2 * (l * w + w * h + h * l)
    total_energy = num_1cm_balls * 1 + num_2cm_balls * 20
    
    # Output the thinking process and the final equation numbers as requested.
    print("Problem Formulation: Yes, this problem can be formulated for an IP solver by using a simplified grid-based packing model.")
    print("\n--- Optimal Solution Search ---")
    print(f"Strategy: Maximize the number of highly efficient 2-cm balls using a box with dimensions L=4a, W=4b, H=4c.")
    print(f"Optimal integer factors found: a={best_abc[0]}, b={best_abc[1]}, c={best_abc[2]}")
    
    print("\n--- Final Configuration ---")
    print(f"Optimal Container: {container_description}")
    print(f"Surface Area Equation: 2 * ({l}*{w} + {w}*{h} + {h}*{l}) = {final_surface_area} cm^2 (<= 1050 cm^2)")
    print(f"Number of 1-cm balls (a): {num_1cm_balls}")
    print(f"Number of 2-cm balls (b): {num_2cm_balls}")
    print(f"Total Energy Equation: {num_1cm_balls} * 1 MJ + {num_2cm_balls} * 20 MJ = {total_energy} MJ")

    # Output the final answer in the specified format.
    final_answer_string = f"[{container_description}]{num_1cm_balls};{num_2cm_balls}"
    print(f"\n<<<Yes>>>{final_answer_string}>>>")

solve_pioneer_probe_packing()