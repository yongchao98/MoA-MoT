import math

def solve_pioneer_probe_packing():
    """
    Solves the exo-planet probe packing problem by searching for the
    optimal container that maximizes the potential energy of packed balls.
    """
    
    max_energy = -1
    best_config = None
    
    # Precision step for container dimensions
    STEP = 0.5
    
    # We will search for the optimal box dimensions (l, w, h).
    # The optimal shape that maximizes volume for a given surface area is a cube.
    # A cube with SA=1050 has sides of sqrt(1050/6) ~= 13.2 cm.
    # Dimensions will be multiples of ball diameters (4 cm), so we search around this size.
    # Search range from 4.0 cm to 25.0 cm should be sufficient.
    dim_range = [i * STEP for i in range(int(4.0/STEP), int(25.0/STEP) + 1)]

    # Iterate through possible dimensions for a box, assuming l >= w >= h
    # to avoid redundant calculations of rotated boxes.
    for l in dim_range:
        for w in dim_range:
            if w > l:
                continue
            for h in dim_range:
                if h > w:
                    continue

                # Constraint 1: Check if the surface area is within the limit
                surface_area = 2 * (l*w + l*h + w*h)
                if surface_area > 1050:
                    continue
                    
                # Using our grid packing model to calculate the number of balls
                
                # Step 1: Pack as many large balls (radius 2cm, diameter 4cm) as possible
                if l < 4 or w < 4 or h < 4:
                    n2 = 0
                else:
                    nx_large = math.floor(l / 4.0)
                    ny_large = math.floor(w / 4.0)
                    nz_large = math.floor(h / 4.0)
                    n2 = nx_large * ny_large * nz_large

                # Step 2: Pack small balls (radius 1cm, diameter 2cm) in the cubic voids
                #          that exist between the large balls.
                nx_small_voids = max(0, nx_large - 1)
                ny_small_voids = max(0, ny_large - 1)
                nz_small_voids = max(0, nz_large - 1)
                n1 = nx_small_voids * ny_small_voids * nz_small_voids

                # Calculate the total energy for this configuration
                current_energy = 20 * n2 + 1 * n1
                
                # Check if this is the best solution found so far
                if current_energy > max_energy:
                    max_energy = current_energy
                    best_config = {
                        'type': 'box',
                        'dims': (l, w, h),
                        'n1': n1,
                        'n2': n2,
                        'energy': current_energy,
                        'sa': surface_area
                    }

    # Format the final answer as specified
    if best_config:
        dims = best_config['dims']
        # Use str() to handle floats like 16.0 becoming "16"
        desc = f"box {str(dims[0])}x{str(dims[1])}x{str(dims[2])}"
        num_1cm_balls = best_config['n1']
        num_2cm_balls = best_config['n2']
        
        # The prompt asks to formulate this problem. Yes, it can be formulated as a complex MINLP problem.
        # This code finds the solution based on a strong heuristic for this problem formulation.
        final_answer = f"[{desc}]{num_1cm_balls};{num_2cm_balls}"
        print(final_answer)
    else:
        # Fallback if no solution is found
        print("[0]")

solve_pioneer_probe_packing()
<<<[box 16.5x16x8]9;32>>>