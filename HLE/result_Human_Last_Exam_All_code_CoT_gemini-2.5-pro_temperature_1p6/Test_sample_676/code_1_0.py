import math

def solve_container_problem():
    """
    Solves the energy ball container optimization problem.
    """
    # --- Constants ---
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = 4.0 # cm
    ENERGY_PER_BALL = 25  # MJ
    COST_PER_BALL = 1000  # USD
    MATERIAL_COST_PER_CM2 = 200  # USD
    REQUIRED_ENERGY = 1000  # MJ
    PRECISION = 0.5 # cm

    # --- Step 1: Minimum number of balls needed ---
    num_balls_needed = math.ceil(REQUIRED_ENERGY / ENERGY_PER_BALL)

    # --- Data structure to hold the best found solution ---
    best_solution = {
        'cost': float('inf'),
        'type': None,
        'num_balls': 0,
        'surface_area': 0
    }

    # Helper function to find integer triplets for a number
    def find_factors(n):
        factors = set()
        for i in range(1, int(n**0.5) + 1):
            if n % i == 0:
                for j in range(i, int((n / i)**0.5) + 1):
                    if (n / i) % j == 0:
                        k = n // (i * j)
                        if i * j * k == n:
                            factors.add(tuple(sorted((i, j, k))))
        return list(factors)

    # Helper function to find integer pairs for a number
    def find_2d_factors(n):
        factors = []
        for i in range(1, int(n**0.5) + 1):
            if n % i == 0:
                factors.append((i, n // i))
        return factors

    # --- Step 2: Box Optimization ---
    # We check for total balls from 40 up to a reasonable limit (e.g., 80).
    # Cost tends to increase with more balls unless a much more efficient packing is found.
    for n_balls in range(num_balls_needed, 81):
        for nx, ny, nz in find_factors(n_balls):
            L = nx * BALL_DIAMETER
            W = ny * BALL_DIAMETER
            H = nz * BALL_DIAMETER
            
            surface_area = 2 * (L*W + W*H + H*L)
            cost_container = surface_area * MATERIAL_COST_PER_CM2
            cost_balls = n_balls * COST_PER_BALL
            total_cost = cost_container + cost_balls

            if total_cost < best_solution['cost']:
                best_solution.update({
                    'cost': total_cost, 'type': 'box', 
                    'num_balls': n_balls, 'surface_area': surface_area
                })

    # --- Step 3: Cylinder Optimization ---
    # We iterate through the number of layers (nz)
    for nz in range(1, num_balls_needed + 1):
        min_n_layer = math.ceil(num_balls_needed / nz)
        # Check a few more balls per layer in case it allows for a better 2D arrangement
        for n_layer in range(min_n_layer, min_n_layer + 10):
            H = nz * BALL_DIAMETER
            
            # Find the best rectangular arrangement (nx_l, ny_l) for n_layer balls
            # that minimizes the diagonal, and thus the containing cylinder's radius.
            best_diag_sq = float('inf')
            best_nx_l, best_ny_l = 0, 0
            for nx_l, ny_l in find_2d_factors(n_layer):
                diag_sq = nx_l**2 + ny_l**2
                if diag_sq < best_diag_sq:
                    best_diag_sq = diag_sq

            # Calculate Radius from diagonal, enforcing precision
            R_uncorrected = 0.5 * BALL_DIAMETER * math.sqrt(best_diag_sq)
            R = math.ceil(R_uncorrected / PRECISION) * PRECISION
            
            surface_area = 2 * math.pi * R * (H + R)
            n_balls_total = nz * n_layer
            cost_container = surface_area * MATERIAL_COST_PER_CM2
            cost_balls = n_balls_total * COST_PER_BALL
            total_cost = cost_container + cost_balls
            
            if total_cost < best_solution['cost']:
                 best_solution.update({
                    'cost': total_cost, 'type': 'cylinder', 
                    'num_balls': n_balls_total, 'surface_area': surface_area
                })


    # --- Step 4: Print Final Result ---
    if best_solution['cost'] == float('inf'):
        print("Could not find a valid design.")
        # This will result in an answer of 0.
    else:
        final_balls = best_solution['num_balls']
        final_sa = best_solution['surface_area']
        final_cost = best_solution['cost']

        print(f"The optimal design is a {best_solution['type']} container.")
        print("The final equation for the total cost C is:")
        # Final cost C = (Number of balls * Cost per ball) + (Surface Area * Material Cost per cm^2)
        print(f"C = {final_balls} * {COST_PER_BALL} + {final_sa:.2f} * {MATERIAL_COST_PER_CM2} = {final_cost:.2f}")


if __name__ == '__main__':
    solve_container_problem()
