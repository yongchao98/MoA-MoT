import math

def solve_pioneer_probe_design():
    """
    Calculates the minimum cost design for a container of energy balls.
    """

    # --- Step 1: Define Constants ---
    E_NEEDED = 1000.0
    E_PER_BALL = 25.0
    COST_PER_BALL = 1000.0
    COST_PER_AREA = 200.0
    R_BALL = 2.0
    D_BALL = 4.0
    PRECISION = 0.5

    def find_best_box_sa(num_balls):
        """
        Finds the minimum surface area for a box-shaped container that can hold a given number of balls.
        Assumes a simple grid packing.
        """
        min_sa = float('inf')
        
        # To minimize surface area for a given volume, the shape should be as close to a cube as possible.
        # We search for integer dimensions (nl, nw, nh) for the ball grid that minimize surface area
        # for a container that can hold at least num_balls.
        limit1 = int(num_balls**(1/3.0)) + 2
        for nl in range(1, limit1):
            limit2 = int((num_balls / nl)**0.5) + 2
            for nw in range(nl, limit2):
                # Calculate the number of layers needed
                nh = math.ceil(num_balls / (nl * nw))
                
                # Calculate container dimensions based on ball grid
                L = nl * D_BALL
                W = nw * D_BALL
                H = nh * D_BALL
                
                # Calculate surface area
                sa = 2 * (L*W + W*H + H*L)
                
                if sa < min_sa:
                    min_sa = sa
                    
        return min_sa

    def find_best_cylinder_sa(num_balls):
        """
        Finds the minimum surface area for a cylindrical container that can hold a given number of balls.
        Assumes balls are stacked in layers, with each layer being a grid.
        """
        min_sa = float('inf')
        
        # Iterate through possible number of layers (nh)
        for nh in range(1, num_balls + 1):
            n_per_layer = math.ceil(num_balls / nh)
            
            # For each nh, find the best (nx, ny) grid arrangement in the base to minimize radius
            min_r_req = float('inf')

            limit_nx = int(n_per_layer**0.5) + 2
            for nx in range(1, limit_nx):
                ny = math.ceil(n_per_layer / nx)
                
                # Required radius to circumscribe the (nx, ny) grid of balls
                r_req = R_BALL * math.sqrt(nx**2 + ny**2)
                if r_req < min_r_req:
                    min_r_req = r_req

            # Actual radius must be a multiple of PRECISION
            R_cyl = math.ceil(min_r_req / PRECISION) * PRECISION
            H_cyl = nh * D_BALL
            
            sa = 2 * math.pi * R_cyl * (R_cyl + H_cyl)
            
            if sa < min_sa:
                min_sa = sa
                
        return min_sa

    # --- Step 2: Determine minimum number of balls ---
    min_balls_needed = math.ceil(E_NEEDED / E_PER_BALL)
    
    min_total_cost = float('inf')
    best_config = {}

    # --- Step 3 & 4: Iterate and find the optimal design ---
    # We search a range of ball numbers. Adding balls increases cost significantly ($1000/ball),
    # so the optimal number is likely very close to the minimum required. A small search
    # range is sufficient to check if a more efficient packing outweighs the cost of extra balls.
    for num_balls in range(int(min_balls_needed), int(min_balls_needed) + 20):
        cost_balls = num_balls * COST_PER_BALL
        
        # --- Box Calculation ---
        sa_box = find_best_box_sa(num_balls)
        cost_box = cost_balls + sa_box * COST_PER_AREA
        
        if cost_box < min_total_cost:
            min_total_cost = cost_box
            best_config = {
                'cost': cost_box,
                'balls': num_balls,
                'sa': sa_box
            }
            
        # --- Cylinder Calculation ---
        sa_cyl = find_best_cylinder_sa(num_balls)
        cost_cyl = cost_balls + sa_cyl * COST_PER_AREA
        
        if cost_cyl < min_total_cost:
            min_total_cost = cost_cyl
            best_config = {
                'cost': cost_cyl,
                'balls': num_balls,
                'sa': sa_cyl
            }

    # --- Step 5: Output the final result ---
    if not best_config:
        print("0")
        print("<<<0>>>")
        return

    final_balls = float(best_config['balls'])
    final_sa = float(best_config['sa'])
    final_total_cost = float(best_config['cost'])
    
    # Print the equation for the total cost calculation
    print(f"{final_balls} * {COST_PER_BALL} + {final_sa} * {COST_PER_AREA} = {final_total_cost}")
    
    # Print the final answer in the required format
    print(f"<<<{final_total_cost}>>>")

solve_pioneer_probe_design()