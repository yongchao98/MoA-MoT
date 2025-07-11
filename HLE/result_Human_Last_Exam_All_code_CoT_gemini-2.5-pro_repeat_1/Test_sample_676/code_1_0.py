import math
import itertools

def solve_container_problem():
    """
    This function solves the container design optimization problem.
    It systematically searches for the configuration (number of balls, arrangement)
    and container type (box or cylinder) that yields the minimum total cost,
    adhering to all specified constraints.
    """
    # --- Problem Constants ---
    BALL_ENERGY = 25  # MJ
    MIN_ENERGY = 1000  # MJ
    BALL_COST = 1000  # usd
    MATERIAL_COST = 200  # usd/cm2
    BALL_DIAMETER = 4.0  # cm
    PRECISION = 0.5  # cm

    MIN_BALLS = math.ceil(MIN_ENERGY / BALL_ENERGY)

    min_total_cost = float('inf')
    best_config = {}

    # --- Search for the optimal configuration ---
    # We search for the number of balls N, starting from the minimum required.
    # The ball cost increases linearly with N, so the optimal N is likely not far from the minimum.
    # We set a search limit for N to 100, which is a reasonable heuristic.
    for n_balls in range(MIN_BALLS, 100):
        # Find all unique sets of 3 integer factors (nx, ny, nz) for n_balls
        factors_set = set()
        for i in range(1, int(n_balls**(1/3.0)) + 2):
            if n_balls % i == 0:
                rem1 = n_balls // i
                for j in range(i, int(rem1**0.5) + 2):
                    if rem1 % j == 0:
                        k = rem1 // j
                        # Check if we found a valid 3-factor combination
                        if i * j * k == n_balls:
                            factors_set.add(tuple(sorted((i, j, k))))

        for factors in factors_set:
            nx, ny, nz = factors
            
            cost_balls = n_balls * BALL_COST

            # The block of spheres has these dimensions.
            # Consider all permutations for the container's L, W, H.
            dims_block = (nx * BALL_DIAMETER, ny * BALL_DIAMETER, nz * BALL_DIAMETER)
            for L, W, H in set(itertools.permutations(dims_block)):

                # --- Case 1: Box Container ---
                sa_box = 2 * (L * W + L * H + W * H)
                total_cost_box = cost_balls + sa_box * MATERIAL_COST

                if total_cost_box < min_total_cost:
                    min_total_cost = total_cost_box
                    best_config = {
                        'cost': total_cost_box,
                        'type': 'Box',
                        'N': n_balls,
                        'arrangement': tuple(sorted((int(L/BALL_DIAMETER), int(W/BALL_DIAMETER), int(H/BALL_DIAMETER)))),
                        'dims': tuple(sorted((L, W, H))),
                        'SA': sa_box,
                    }

                # --- Case 2: Cylinder Container ---
                # The base of the cylinder must enclose a rectangle of L x W
                radius_raw = 0.5 * math.sqrt(L**2 + W**2)
                # Apply the manufacturing precision constraint
                radius_cyl = math.ceil(radius_raw / PRECISION) * PRECISION
                height_cyl = H

                sa_cyl = 2 * math.pi * radius_cyl**2 + 2 * math.pi * radius_cyl * height_cyl
                total_cost_cyl = cost_balls + sa_cyl * MATERIAL_COST

                if total_cost_cyl < min_total_cost:
                    min_total_cost = total_cost_cyl
                    best_config = {
                        'cost': total_cost_cyl,
                        'type': 'Cylinder',
                        'N': n_balls,
                        'arrangement': tuple(sorted((int(L/BALL_DIAMETER), int(W/BALL_DIAMETER), int(H/BALL_DIAMETER)))),
                        'base_dims': tuple(sorted((L, W))),
                        'dims': (radius_cyl, height_cyl),
                        'SA': sa_cyl,
                    }

    # --- Print the final results of the best found design ---
    if not best_config:
        print("No solution found.")
        return

    N = best_config['N']
    cost_balls = N * BALL_COST
    cost_container = best_config['cost'] - cost_balls

    print("This script calculates the lowest total cost for a container and energy balls.")
    print("It considers both box and cylinder shapes, optimizing for the number of balls and container dimensions.")
    print("\n--- Optimal Design Found ---")
    
    if best_config['type'] == 'Box':
        L, W, H = best_config['dims']
        nx, ny, nz = best_config['arrangement']
        print(f"Container Type: Box")
        print(f"Number of Balls (N): {N}")
        print(f"Ball Arrangement: {nx} x {ny} x {nz} grid")
        print(f"\nFinal Equation:")
        print(f"Cost_balls = {N} * {BALL_COST} = {cost_balls:.0f}")
        print(f"Container_dimensions = ({L:.1f}, {W:.1f}, {H:.1f}) cm")
        print(f"Container_SA = 2 * ({L:.1f}*{W:.1f} + {L:.1f}*{H:.1f} + {W:.1f}*{H:.1f}) = {best_config['SA']:.2f} cm2")
        print(f"Cost_container = {best_config['SA']:.2f} * {MATERIAL_COST} = {cost_container:.0f}")
        print(f"Total_Cost = {cost_balls:.0f} + {cost_container:.0f} = {best_config['cost']:.0f}")

    elif best_config['type'] == 'Cylinder':
        R, H = best_config['dims']
        base_L, base_W = best_config['base_dims']
        print(f"Container Type: Cylinder")
        print(f"Number of Balls (N): {N}")
        print(f"Ball Arrangement: Layers of spheres on a {base_L:.1f}x{base_W:.1f} cm base, stacked to a height of {H:.1f} cm")
        print(f"\nFinal Equation:")
        print(f"Cost_balls = {N} * {BALL_COST} = {cost_balls:.0f}")
        print(f"Container_dimensions = (Radius={R:.1f}, Height={H:.1f}) cm")
        print(f"Container_SA = 2*pi*{R:.1f}^2 + 2*pi*{R:.1f}*{H:.1f} = {best_config['SA']:.2f} cm2")
        print(f"Cost_container = {best_config['SA']:.2f} * {MATERIAL_COST} = {cost_container:.0f}")
        print(f"Total_Cost = {cost_balls:.0f} + {cost_container:.0f} = {best_config['cost']:.0f}")

solve_container_problem()