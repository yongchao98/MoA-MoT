import math

def solve_pioneer_probe_design():
    """
    Finds the optimal container design (box or cylinder) to pack energy balls
    with the lowest total cost, subject to energy and material constraints.
    """

    # --- Constants ---
    R_BALL = 2.0  # cm
    D_BALL = 4.0  # cm
    ENERGY_PER_BALL = 30  # MJ
    COST_PER_BALL = 1000  # USD
    SA_MAX = 1000.0  # cm^2
    COST_PER_CM2 = 200.0  # USD
    E_TOTAL_MIN = 1000.0  # MJ
    PRECISION = 0.5  # cm

    # --- Derived Requirements ---
    MIN_BALLS_NEEDED = math.ceil(E_TOTAL_MIN / ENERGY_PER_BALL)
    D_BALL_UNITS = int(D_BALL / PRECISION)

    min_total_cost = float('inf')
    best_design = None

    # --- Part 1: Search for the best Box Container ---
    # Search range for dimensions in PRECISION units.
    # A cube with SA=1000 has a side of ~12.9cm -> 26 units.
    # We search slightly higher to be safe.
    max_dim_unit_box = 35 

    for l_unit in range(1, max_dim_unit_box + 1):
        for w_unit in range(1, l_unit + 1):
            for h_unit in range(1, w_unit + 1):
                sa = (l_unit * w_unit + w_unit * h_unit + h_unit * l_unit) * 2 * (PRECISION**2)
                
                if sa > SA_MAX:
                    if h_unit == 1 and w_unit == 1: # No point increasing l_unit if this fails
                         break 
                    if h_unit == 1: # No point increasing w_unit if this fails
                         break
                    continue

                nx = l_unit // D_BALL_UNITS
                ny = w_unit // D_BALL_UNITS
                nz = h_unit // D_BALL_UNITS
                num_balls = nx * ny * nz

                if num_balls < MIN_BALLS_NEEDED:
                    continue

                cost_container = sa * COST_PER_CM2
                cost_balls = num_balls * COST_PER_BALL
                total_cost = cost_container + cost_balls

                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        'type': 'Box',
                        'dims': (l_unit * PRECISION, w_unit * PRECISION, h_unit * PRECISION),
                        'num_balls': num_balls,
                        'sa': sa,
                        'total_cost': total_cost,
                    }

    # --- Part 2: Search for the best Cylinder Container ---
    # Max radius for SA<=1000 for a flat cylinder (H=0) is R~12.6cm -> 25 units.
    max_r_unit_cyl = 26

    for r_unit in range(1, max_r_unit_cyl + 1):
        R = r_unit * PRECISION
        if R < R_BALL:
            continue

        # Calculate balls per layer using a grid packing model
        n_per_layer = 0
        if R >= R_BALL:
            max_i = int((R - R_BALL) / D_BALL)
            for i in range(-max_i, max_i + 1):
                for j in range(-max_i, max_i + 1):
                    x, y = i * D_BALL, j * D_BALL
                    if math.sqrt(x**2 + y**2) <= R - R_BALL:
                        n_per_layer += 1
        
        if n_per_layer == 0:
            continue
        
        # Search for height h_unit
        # From 2*pi*R*(R+H) <= SA_MAX => H <= SA_MAX/(2*pi*R) - R
        max_h_cm = (SA_MAX / (2 * math.pi * R)) - R if R > 0 else 0
        if max_h_cm <= 0: continue
        max_h_unit = int(max_h_cm / PRECISION)

        for h_unit in range(1, max_h_unit + 1):
            H = h_unit * PRECISION
            sa = 2 * math.pi * R * (R + H)

            n_layers = int(H // D_BALL)
            num_balls = n_layers * n_per_layer

            if num_balls < MIN_BALLS_NEEDED:
                continue

            cost_container = sa * COST_PER_CM2
            cost_balls = num_balls * COST_PER_BALL
            total_cost = cost_container + cost_balls
            
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design = {
                    'type': 'Cylinder',
                    'dims': (R, H),
                    'num_balls': num_balls,
                    'sa': sa,
                    'total_cost': total_cost,
                }

    # --- Part 3: Output the Final Result ---
    if min_total_cost == float('inf') or best_design is None:
        print(0)
        print("<<<0>>>")
    else:
        container_type = best_design['type']
        dims = best_design['dims']
        num_balls = best_design['num_balls']
        sa = best_design['sa']
        
        cost_container = sa * COST_PER_CM2
        cost_balls = num_balls * COST_PER_BALL
        
        print(f"Optimal design is a {container_type} container.")
        if container_type == 'Box':
            print(f"Dimensions (L x W x H): {dims[0]} cm x {dims[1]} cm x {dims[2]} cm")
        else:
            print(f"Dimensions (Radius x Height): {dims[0]} cm x {dims[1]} cm")
        
        print(f"The container will hold {num_balls} energy balls, providing {num_balls * ENERGY_PER_BALL} MJ of energy.")
        print(f"The surface area is {sa:.2f} cm^2 and the total cost is ${min_total_cost:.2f}.")
        print("\nCost breakdown:")
        print(f"Final equation: {cost_container:.2f} ({sa:.2f} * {COST_PER_CM2}) + {cost_balls:.2f} ({num_balls} * {COST_PER_BALL}) = {min_total_cost:.2f}")

        final_answer = int(round(min_total_cost))
        print(f"<<<{final_answer}>>>")

solve_pioneer_probe_design()