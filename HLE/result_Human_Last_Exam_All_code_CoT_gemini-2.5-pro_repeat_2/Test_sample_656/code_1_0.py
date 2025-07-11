import math

def solve_packing_problem():
    """
    Finds the optimal box container to maximize energy from packed spheres.
    """
    max_surface_area = 1050.0
    
    # Ball properties
    # We only consider 2-cm balls as they are far more energy-dense.
    ball_radius = 2.0
    ball_diameter = 2 * ball_radius
    energy_per_ball = 20

    # Search parameters
    dim_step = 0.5
    # Max possible dimension for a cube is sqrt(1050/6) ~= 13.2 cm.
    # We search a bit higher for non-cubic shapes.
    # Dimensions are i * dim_step, so max_i of 60 corresponds to 30 cm.
    max_i = 61 

    best_energy = 0
    best_config = {}

    # Iterate through possible dimensions l, w, h as multiples of 0.5 cm
    # We use i, j, k as integer multiples of dim_step
    # To avoid permutations, we enforce i <= j <= k
    for i in range(1, max_i):
        l = i * dim_step
        if l < ball_diameter:
            continue

        for j in range(i, max_i):
            w = j * dim_step
            if w < ball_diameter:
                continue

            # Optimization: if the area of the base alone is too large, stop.
            if 2 * l * w > max_surface_area:
                break
            
            for k in range(j, max_i):
                h = k * dim_step
                if h < ball_diameter:
                    continue

                # Calculate surface area
                surface_area = 2 * (l * w + l * h + w * h)

                # If SA exceeds the limit, any larger h will also exceed it.
                if surface_area > max_surface_area:
                    break
                
                # Calculate number of balls using simple lattice packing
                num_balls_l = math.floor(l / ball_diameter)
                num_balls_w = math.floor(w / ball_diameter)
                num_balls_h = math.floor(h / ball_diameter)
                
                total_balls = num_balls_l * num_balls_w * num_balls_h
                current_energy = total_balls * energy_per_ball

                if current_energy > best_energy:
                    best_energy = current_energy
                    best_config = {
                        "shape": "box",
                        "l": l,
                        "w": w,
                        "h": h,
                        "n1": 0,
                        "n2": total_balls,
                        "surface_area": surface_area
                    }
                # Tie-breaker: if energies are equal, prefer the one with larger SA
                # This means it uses the material more efficiently for the same result.
                elif current_energy == best_energy and current_energy > 0:
                    if surface_area > best_config.get("surface_area", 0):
                         best_config = {
                            "shape": "box",
                            "l": l,
                            "w": w,
                            "h": h,
                            "n1": 0,
                            "n2": total_balls,
                            "surface_area": surface_area
                        }

    # Format the final output string
    if best_config:
        l, w, h = best_config["l"], best_config["w"], best_config["h"]
        # Helper to format numbers to remove trailing .0
        def fmt(n):
            return int(n) if n == int(n) else n
        
        container_desc = f"box {fmt(l)}x{fmt(w)}x{fmt(h)}"
        num_1cm_balls = best_config["n1"]
        num_2cm_balls = best_config["n2"]
        
        result_string = f"[{container_desc}]{num_1cm_balls};{num_2cm_balls}"
        print(result_string)
    else:
        print("[0]")

solve_packing_problem()
