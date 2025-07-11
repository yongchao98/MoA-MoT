import math

def solve_container_problem():
    """
    Solves the container optimization problem by finding the lowest cost
    design for either a box or a cylinder container.
    """
    # --- Problem Constants ---
    REQUIRED_ENERGY_MJ = 1000
    ENERGY_PER_BALL_MJ = 25
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    BALL_COST_USD = 1000
    MATERIAL_COST_USD_PER_CM2 = 200
    PRECISION_CM = 0.5

    # --- Basic Calculations ---
    min_balls_needed = math.ceil(REQUIRED_ENERGY_MJ / ENERGY_PER_BALL_MJ)

    # --- Box Optimization ---
    best_box = {'cost': float('inf')}

    def find_factors_triplets(n):
        """Finds all unique integer triplets (i, j, k) such that i*j*k = n."""
        factors = []
        for i in range(1, int(n**(1/3.0)) + 2):
            if n % i == 0:
                for j in range(i, int(math.sqrt(n / i)) + 2):
                    if (n / i) % j == 0:
                        k = n // (i * j)
                        if i * j * k == n:
                            factors.append(tuple(sorted((i, j, k))))
        return sorted(list(set(factors)))

    # We test a range of ball counts starting from the minimum required.
    for num_balls in range(min_balls_needed, min_balls_needed + 21):
        triplets = find_factors_triplets(num_balls)
        if not triplets:  # Prime number case
             triplets.append((1,1,num_balls))

        for nx, ny, nz in triplets:
            L = nx * BALL_DIAMETER_CM
            W = ny * BALL_DIAMETER_CM
            H = nz * BALL_DIAMETER_CM

            surface_area = 2 * (L*W + L*H + W*H)
            material_cost = surface_area * MATERIAL_COST_USD_PER_CM2
            ball_cost = num_balls * BALL_COST_USD
            total_cost = material_cost + ball_cost

            if total_cost < best_box['cost']:
                best_box = {
                    'cost': total_cost, 'type': 'Box', 'num_balls': num_balls,
                    'arrangement': (nx, ny, nz), 'dims': (L, W, H),
                    'surface_area': surface_area, 'material_cost': material_cost,
                    'ball_cost': ball_cost
                }

    # --- Cylinder Optimization ---
    # R_unit(k) is the radius of the smallest circle that can contain k unit circles (radius 1).
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cci/cci.html
    R_unit_k = {
        1: 1.000, 2: 2.000, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.000, 7: 3.000,
        8: 3.304, 9: 3.613, 10: 3.813, 11: 3.924, 12: 4.030, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 4.864,
        20: 5.122, 21: 5.248, 22: 5.341, 23: 5.534, 24: 5.607, 25: 5.760,
        26: 5.869, 27: 5.869, 28: 6.091, 29: 6.166, 30: 6.273, 31: 6.402,
        32: 6.467, 33: 6.593, 34: 6.666, 35: 6.748, 36: 6.864, 37: 6.864,
        38: 7.070, 39: 7.141, 40: 7.228
    }

    best_cylinder = {'cost': float('inf')}

    for nz in range(1, min_balls_needed + 1):
        k_per_layer = math.ceil(min_balls_needed / nz)
        num_balls = nz * k_per_layer

        H = nz * BALL_DIAMETER_CM
        
        min_radius_for_packing = R_unit_k.get(k_per_layer, math.sqrt(k_per_layer / 0.9069)) * BALL_RADIUS_CM
        R = math.ceil(min_radius_for_packing / PRECISION_CM) * PRECISION_CM
        
        surface_area = 2 * math.pi * R**2 + 2 * math.pi * R * H
        material_cost = surface_area * MATERIAL_COST_USD_PER_CM2
        ball_cost = num_balls * BALL_COST_USD
        total_cost = material_cost + ball_cost
        
        if total_cost < best_cylinder['cost']:
            best_cylinder = {
                'cost': total_cost, 'type': 'Cylinder', 'num_balls': num_balls,
                'arrangement': (nz, k_per_layer), 'dims': (R, H),
                'surface_area': surface_area, 'material_cost': material_cost,
                'ball_cost': ball_cost
            }

    # --- Comparison and Result ---
    winner = best_box if best_box['cost'] < best_cylinder['cost'] else best_cylinder
        
    final_cost = winner['cost']
    
    print("--- Optimal Container Design ---")
    print(f"Container Type: {winner['type']}")
    if winner['type'] == 'Box':
        nx, ny, nz = winner['arrangement']
        L, W, H = winner['dims']
        print(f"Ball Arrangement: {nx} x {ny} x {nz}")
        print(f"Container Dimensions (cm): Length={L}, Width={W}, Height={H}")
    else:  # Cylinder
        nz, k = winner['arrangement']
        R, H = winner['dims']
        print(f"Ball Arrangement: {nz} layers of {k} balls each")
        print(f"Container Dimensions (cm): Radius={R}, Height={H}")

    print(f"\nNumber of energy balls: {winner['num_balls']}")
    print(f"Surface Area: {winner['surface_area']:.2f} cm^2")
    
    print("\n--- Cost Calculation ---")
    # Print the numbers in the final equation as requested
    print("Cost of Energy Balls:")
    print(f"{winner['num_balls']} balls * ${BALL_COST_USD}/ball = ${winner['ball_cost']:.2f}")
    
    print("\nCost of Container Material:")
    print(f"{winner['surface_area']:.2f} cm^2 * ${MATERIAL_COST_USD_PER_CM2}/cm^2 = ${winner['material_cost']:.2f}")
    
    print("\n--- Final Total Cost Equation ---")
    print(f"Total Cost C = (Number of Balls * Cost per Ball) + (Surface Area * Cost per cm^2)")
    print(f"C = ({winner['num_balls']} * {BALL_COST_USD}) + ({winner['surface_area']:.2f} * {MATERIAL_COST_USD_PER_CM2})")
    print(f"C = ${winner['ball_cost']:.2f} + ${winner['material_cost']:.2f}")
    print(f"C = ${winner['cost']:.2f}")

    print(f"\nConclusion: The lowest total cost is for the {winner['type']} container.")
    
    print(f"<<<{final_cost:.0f}>>>")

# Execute the function to find and print the solution
solve_container_problem()
