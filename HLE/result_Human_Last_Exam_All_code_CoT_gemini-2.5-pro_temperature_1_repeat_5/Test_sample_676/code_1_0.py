import math

def solve_container_problem():
    """
    Solves the container optimization problem to find the minimum cost design.
    """
    # 1. Define Problem Constants
    ENERGY_PER_BALL = 25  # MJ
    MIN_TOTAL_ENERGY = 1000  # MJ
    BALL_COST = 1000  # usd
    BALL_RADIUS = 2  # cm
    BALL_DIAMETER = 4 # cm
    CONTAINER_MATERIAL_COST = 200  # usd per cm2
    PRECISION = 0.5  # cm

    # Data for circle packing (R_container / R_ball)
    # Source: "Packing of equal circles in a circle" by E. Specht and others.
    PACKING_RATIOS = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 4.864,
        20: 5.122, 21: 5.248, 22: 5.340, 23: 5.515, 24: 5.606, 25: 5.680,
        26: 5.865, 27: 5.889, 28: 6.000, 29: 6.102, 30: 6.126, 31: 6.264,
        32: 6.366, 33: 6.442, 34: 6.505, 35: 6.634, 36: 6.666, 37: 6.747,
        38: 6.848, 39: 6.942, 40: 7.013, 41: 7.143, 42: 7.218, 43: 7.318,
        44: 7.375, 45: 7.456, 46: 7.505, 47: 7.625, 48: 7.653, 49: 7.712,
        50: 7.785
    }

    min_balls_required = math.ceil(MIN_TOTAL_ENERGY / ENERGY_PER_BALL)
    min_total_cost = float('inf')
    best_design = {}

    # 2. Optimize the Box Container
    # Search for number of balls n from min_balls_required up to a reasonable limit
    # The cost of balls increases linearly, so the optimal n should be close to min_balls_required
    search_limit_n = min_balls_required + 25 
    for n in range(min_balls_required, search_limit_n):
        limit1 = int(n**(1/3.0)) + 2
        for nx in range(1, limit1):
            if n % nx == 0:
                n_rem = n // nx
                limit2 = int(math.sqrt(n_rem)) + 2
                for ny in range(nx, limit2):
                    if n_rem % ny == 0:
                        nz = n_rem // ny
                        if nz >= ny: # Ensures nx <= ny <= nz to avoid permutations
                            cost_balls = n * BALL_COST
                            L, W, H = nx * BALL_DIAMETER, ny * BALL_DIAMETER, nz * BALL_DIAMETER
                            surface_area = 2 * (L*W + L*H + W*H)
                            cost_container = surface_area * CONTAINER_MATERIAL_COST
                            total_cost = cost_balls + cost_container
                            
                            if total_cost < min_total_cost:
                                min_total_cost = total_cost
                                best_design = {
                                    'type': 'Box', 'cost': total_cost, 'n_balls': n,
                                    'arrangement': (nx, ny, nz), 'dimensions': (L, W, H),
                                    'cost_balls': cost_balls, 'cost_container': cost_container
                                }

    # 3. Optimize the Cylinder Container
    search_limit_k = min_balls_required + 11 
    for k in range(1, search_limit_k):
        if k not in PACKING_RATIOS:
            continue
            
        nh = math.ceil(min_balls_required / k)
        n = nh * k
        
        cost_balls = n * BALL_COST
        H = nh * BALL_DIAMETER
        min_packing_radius = PACKING_RATIOS[k] * BALL_RADIUS
        R = math.ceil(min_packing_radius / PRECISION) * PRECISION
        
        surface_area = 2 * math.pi * R**2 + 2 * math.pi * R * H
        cost_container = surface_area * CONTAINER_MATERIAL_COST
        total_cost = cost_balls + cost_container
        
        if total_cost < min_total_cost:
            min_total_cost = total_cost
            best_design = {
                'type': 'Cylinder', 'cost': total_cost, 'n_balls': n,
                'arrangement': (k, nh), 'dimensions': (R, H),
                'cost_balls': cost_balls, 'cost_container': cost_container
            }

    # 4. Compare and Finalize
    print("This task is to design a container to pack energy balls for a total of at least 1000 MJ with minimum total cost.")
    print(f"Each ball provides {ENERGY_PER_BALL} MJ, so we need at least {min_balls_required} balls.")
    print("Total Cost = (Number of Balls * $1,000) + (Container Surface Area cm^2 * $200).")
    print("-" * 30)
    
    if not best_design:
        final_cost = 0
        print("No solution found.")
    else:
        final_cost = best_design['cost']
        print(f"The optimal design found is a {best_design['type']} container.")
        
        if best_design['type'] == 'Box':
            nx, ny, nz = best_design['arrangement']
            L, W, H = best_design['dimensions']
            print(f"It packs {best_design['n_balls']} balls in a {nx}x{ny}x{nz} arrangement.")
            print(f"The container dimensions are L={L} cm, W={W} cm, H={H} cm.")
            print("\nFinal Cost Calculation:")
            print(f"Total Cost = ({best_design['n_balls']} * {BALL_COST}) + (2 * ({L}*{W} + {L}*{H} + {W}*{H})) * {CONTAINER_MATERIAL_COST}")
            print(f"Total Cost = ${best_design['cost_balls']} + ${best_design['cost_container']:.2f} = ${final_cost:.2f}")

        elif best_design['type'] == 'Cylinder':
            k, nh = best_design['arrangement']
            R, H = best_design['dimensions']
            print(f"It packs {best_design['n_balls']} balls in {nh} layers of {k} balls each.")
            print(f"The container dimensions are Radius={R} cm, Height={H} cm.")
            print("\nFinal Cost Calculation:")
            print(f"Total Cost = ({best_design['n_balls']} * {BALL_COST}) + (2 * 3.14159 * {R:.1f}**2 + 2 * 3.14159 * {R:.1f} * {H}) * {CONTAINER_MATERIAL_COST}")
            print(f"Total Cost = ${best_design['cost_balls']} + ${best_design['cost_container']:.2f} = ${final_cost:.2f}")

    # Output the final answer in the required format
    print(f"\nThe lowest total cost C is {round(final_cost)}.")
    print(f"<<<{round(final_cost)}>>>")

solve_container_problem()