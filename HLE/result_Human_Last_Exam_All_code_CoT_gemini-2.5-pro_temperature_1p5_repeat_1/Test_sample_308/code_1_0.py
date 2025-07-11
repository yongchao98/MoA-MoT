import math
import itertools

def solve_rocket_problem():
    """
    Finds the cheapest rocket configuration for a given space mission.
    """
    # Step 1: Define component data and calculate required delta-v
    motors = {
        1: {'price': 8000, 'mass': 1.3, 'v_e': 1200},
        2: {'price': 16000, 'mass': 1.54, 'v_e': 2300}
    }

    tanks = {
        1: {'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        2: {'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        3: {'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        4: {'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    }

    # Delta-v required = transfer (Earth->Helioxis) + landing (Helioxis orbit->surface)
    required_dv = 271 + 278

    # Step 2: Generate all possible rocket configurations
    tank_ids = list(tanks.keys())
    three_tank_combos = list(itertools.combinations(tank_ids, 3))
    four_tank_combos = list(itertools.combinations(tank_ids, 4))
    all_tank_combos = three_tank_combos + four_tank_combos

    capable_rockets = []

    # Step 3 & 4: Evaluate each configuration and filter capable ones
    for motor_id, motor_data in motors.items():
        for combo in all_tank_combos:
            # Calculate properties for the current configuration
            total_cost = motor_data['price'] + sum(tanks[tid]['price'] for tid in combo)
            m0 = motor_data['mass'] + sum(tanks[tid]['wet_mass'] for tid in combo)
            mf = motor_data['mass'] + sum(tanks[tid]['dry_mass'] for tid in combo)
            v_e = motor_data['v_e']

            # Calculate achievable delta-v
            if m0 > mf:
                achieved_dv = v_e * math.log(m0 / mf)
            else:
                achieved_dv = 0

            # Check if it meets the requirement
            if achieved_dv >= required_dv:
                capable_rockets.append({
                    'motor': motor_id,
                    'tanks': combo,
                    'cost': total_cost,
                    'v_e': v_e,
                    'm0': m0,
                    'mf': mf,
                    'dv': achieved_dv
                })

    # Step 5: Find the cheapest solution
    if not capable_rockets:
        print("No capable rocket configuration found.")
        return

    cheapest_rocket = min(capable_rockets, key=lambda x: x['cost'])

    # Step 6: Format and print the output
    motor_str = f"({cheapest_rocket['motor']})"
    tanks_str = ", ".join(map(str, sorted(cheapest_rocket['tanks'])))
    answer_format = f"{motor_str} {tanks_str}"

    print(f"The cheapest rocket configuration is: {answer_format}")
    print("\nCalculation for this rocket using the Tsiolkovsky rocket equation:")
    print(f"delta_v = exhaust_velocity * ln(wet_mass / dry_mass)")
    print(f"delta_v = {cheapest_rocket['v_e']} * ln({cheapest_rocket['m0']:.2f} / {cheapest_rocket['mf']:.2f})")
    print(f"Resulting delta_v: {cheapest_rocket['dv']:.2f} m/s (Required: {required_dv} m/s)")


solve_rocket_problem()