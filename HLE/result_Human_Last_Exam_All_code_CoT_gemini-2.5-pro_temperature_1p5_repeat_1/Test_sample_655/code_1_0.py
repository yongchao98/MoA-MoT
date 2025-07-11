import math
import numpy as np

def solve():
    """
    Finds the optimal box container and ball packing to maximize energy.
    """
    max_surface_area = 1050
    best_energy = 0
    best_config = {}

    # Iterate through possible dimensions (l, w, h) which are multiples of 0.5
    # To reduce search space, assume l <= w <= h
    # A single dimension can't be much larger than sqrt(1050/2) ~= 22.9
    # We set a practical limit for the search, e.g., 30 cm
    limit = int(30 / 0.5) 
    
    for i in range(1, limit + 1):
        l = i * 0.5
        for j in range(i, limit + 1):
            w = j * 0.5
            
            # Estimate max possible h to prune search
            # If l*w is large, h must be small. 2*l*w <= max_surface_area
            if 2 * l * w > max_surface_area:
                break
                
            # derived from sa = 2*(lw+lh+wh) => h = (sa/2 - lw) / (l+w)
            h_max = (max_surface_area / 2.0 - l * w) / (l + w)
            if h_max < w: # since we search h >= w
                continue
            
            h_limit = int(h_max / 0.5)

            for k in range(j, h_limit + 1):
                h = k * 0.5

                surface_area = 2 * (l * w + l * h + w * h)

                if surface_area <= max_surface_area:
                    # Calculate the number of balls using the combined lattice packing strategy
                    # First, calculate max number of big balls (2-cm radius)
                    num_b_balls = math.floor(l / 4) * math.floor(w / 4) * math.floor(h / 4)
                    
                    # Then calculate the total number of small ball (1-cm radius) slots
                    total_a_slots = math.floor(l / 2) * math.floor(w / 2) * math.floor(h / 2)

                    # A 2-cm ball takes the space of 8 1-cm balls (2x2x2)
                    num_a_balls = total_a_slots - num_b_balls * 8
                    
                    # Total energy is E = 10*b + 1*a
                    total_energy = num_b_balls * 10 + num_a_balls * 1

                    if total_energy > best_energy:
                        best_energy = total_energy
                        best_config = {
                            "l": l, "w": w, "h": h,
                            "a": num_a_balls, "b": num_b_balls,
                            "sa": surface_area
                        }

    # Print the result
    l, w, h = best_config["l"], best_config["w"], best_config["h"]
    a, b = best_config["a"], best_config["b"]
    energy = best_energy
    
    print(f"Optimal Container: box {l}x{w}x{h}")
    print(f"Surface Area: {best_config['sa']:.2f} cm^2")
    print(f"Number of 1-cm balls (a): {a}")
    print(f"Number of 2-cm balls (b): {b}")
    print(f"Final Equation: {a} * 1 + {b} * 10 = {energy} MJ")
    print(f"<<<[box {l}x{w}x{h}]{a};{b}>>>")

solve()