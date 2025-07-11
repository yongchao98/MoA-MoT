import math

def solve_pioneer_probe_packing():
    """
    Solves the energy ball packing problem for the Pioneer probe.
    This function implements a greedy packing algorithm to find the optimal
    container and ball configuration to maximize stored energy.
    """

    # --- Helper Functions for Geometry and Checks ---

    def check_overlap(center1, radius1, placed_balls):
        """Check if a new sphere overlaps with any in a list of placed balls."""
        x1, y1, z1 = center1
        # Use squared distances to avoid costly square root operations
        radius1_sq = radius1**2
        for center2, radius2 in placed_balls:
            x2, y2, z2 = center2
            dist_sq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
            min_dist = radius1 + radius2
            if dist_sq < min_dist**2 - 1e-9: # Add tolerance for float precision
                return True  # Overlap found
        return False  # No overlap

    # --- Main Packer Function ---

    def pack_container(container_type, dims):
        """
        Performs a greedy grid-based packing for a given container.
        Returns (n1, n2): number of 1cm and 2cm radius balls.
        """
        placed_balls = []
        step = 0.5

        # Define containment check function and grid boundaries based on type
        if container_type == 'box':
            l, w, h = dims['l'], dims['w'], dims['h']
            def is_contained(center, radius):
                x, y, z = center
                return (radius <= x <= l - radius and
                        radius <= y <= w - radius and
                        radius <= z <= h - radius)
            x_range = [i * step for i in range(int(l / step) + 1)]
            y_range = [i * step for i in range(int(w / step) + 1)]
            z_range = [i * step for i in range(int(h / step) + 1)]

        elif container_type == 'sphere':
            R = dims['r']
            R_sq = R**2
            def is_contained(center, radius):
                x, y, z = center
                # Center of the ball must be within a sphere of radius (R - radius)
                return (x**2 + y**2 + z**2) <= (R - radius)**2 + 1e-9
            coord_range = [i * step for i in range(-int(R / step), int(R / step) + 1)]
            x_range, y_range, z_range = coord_range, coord_range, coord_range

        elif container_type == 'cylinder':
            r_cyl, h_cyl = dims['r'], dims['h']
            def is_contained(center, radius):
                x, y, z = center
                # Check radial and height containment
                return ((x**2 + y**2) <= (r_cyl - radius)**2 + 1e-9 and
                        radius <= z <= h_cyl - radius)
            radial_range = [i * step for i in range(-int(r_cyl / step), int(r_cyl / step) + 1)]
            x_range, y_range = radial_range, radial_range
            z_range = [i * step for i in range(int(h_cyl / step) + 1)]
        else:
            return 0, 0

        # --- Packing Simulation ---
        # Phase 1: Pack 2-cm balls (radius=2.0)
        ball_radius = 2.0
        for z in z_range:
            for y in y_range:
                for x in x_range:
                    center = (x, y, z)
                    if is_contained(center, ball_radius):
                        if not check_overlap(center, ball_radius, placed_balls):
                            placed_balls.append((center, ball_radius))
        
        n2 = len(placed_balls)

        # Phase 2: Pack 1-cm balls (radius=1.0)
        ball_radius = 1.0
        for z in z_range:
            for y in y_range:
                for x in x_range:
                    center = (x, y, z)
                    if is_contained(center, ball_radius):
                        if not check_overlap(center, ball_radius, placed_balls):
                            placed_balls.append((center, ball_radius))

        n1 = len(placed_balls) - n2
        return n1, n2

    # --- Main Execution Logic ---
    print("Starting optimization. This process may take a few minutes...")
    
    # Candidate containers that respect SA <= 1050 cm^2
    candidates = [
        {'type': 'sphere', 'dims': {'r': 9.0}, 'desc_raw': 'sphere r=9.0', 'sa': 4*math.pi*9**2},
        {'type': 'box', 'dims': {'l': 13.0, 'w': 13.0, 'h': 13.0}, 'desc_raw': 'box 13.0x13.0x13.0', 'sa': 6*13**2},
        {'type': 'box', 'dims': {'l': 12.0, 'w': 12.0, 'h': 15.5}, 'desc_raw': 'box 12.0x12.0x15.5', 'sa': 2*(12*12 + 12*15.5 + 12*15.5)},
        {'type': 'cylinder', 'dims': {'r': 5.0, 'h': 28.0}, 'desc_raw': 'cylinder r=5.0, h=28.0', 'sa': 2*math.pi*5**2 + 2*math.pi*5*28.0}
    ]

    max_energy = -1
    best_config = None

    for cand in candidates:
        print(f"\nTesting container: {cand['desc_raw']} (SA: {cand['sa']:.2f} cm^2)")
        n1, n2 = pack_container(cand['type'], cand['dims'])
        energy = 10 * n2 + n1
        
        print(f" -> Result: {n2} large balls (2-cm), {n1} small balls (1-cm).")
        print(f" -> Total Energy: {energy} MJ")
        
        if energy > max_energy:
            max_energy = energy
            best_config = {
                'type': cand['type'],
                'dims': cand['dims'],
                'n1': n1,
                'n2': n2,
                'energy': energy
            }

    # --- Format and Print the Final Answer ---
    print("\n------------------------------------")
    print("Optimization Complete. Best configuration found:")
    
    final_n1 = best_config['n1']
    final_n2 = best_config['n2']
    
    # Create the formatted description string like "box 10x10x10"
    def format_dim(d):
        return str(int(d)) if d == int(d) else str(d)

    if best_config['type'] == 'box':
        dims = best_config['dims']
        l, w, h = format_dim(dims['l']), format_dim(dims['w']), format_dim(dims['h'])
        formatted_desc = f"box {l}x{w}x{h}"
    elif best_config['type'] == 'sphere':
        dims = best_config['dims']
        r = format_dim(dims['r'])
        formatted_desc = f"sphere r={r}"
    elif best_config['type'] == 'cylinder':
        dims = best_config['dims']
        r, h = format_dim(dims['r']), format_dim(dims['h'])
        formatted_desc = f"cylinder r={r}, h={h}"

    print(f"Container Description: {formatted_desc}")
    print(f"Number of 1-cm balls (a): {final_n1}")
    print(f"Number of 2-cm balls (b): {final_n2}")
    
    print("\nFinal Energy Calculation:")
    print(f"({final_n2} * 10 MJ) + ({final_n1} * 1 MJ) = {best_config['energy']} MJ")

    final_answer_string = f"[{formatted_desc}]{final_n1};{final_n2}"
    print("\nFinal Answer in required format:")
    print(f'<<<{final_answer_string}>>>')

# Execute the solver
solve_pioneer_probe_packing()