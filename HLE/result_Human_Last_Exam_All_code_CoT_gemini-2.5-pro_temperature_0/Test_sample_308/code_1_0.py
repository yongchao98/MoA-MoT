import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis.
    """

    # Data from the tables
    motors = {
        1: {'price': 8000, 'mass': 1.3, 've': 1200},  # ve in m/s
        2: {'price': 16000, 'mass': 1.54, 've': 2300}
    }

    tanks = {
        1: {'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        2: {'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        3: {'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        4: {'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    }

    # Step 1: Calculate required Delta-V for the mission
    # Low Earth orbit to escape Earth: 28 m/s
    # Earth to Helioxis transfer: 271 m/s
    # Helioxis orbit insertion (same as surface to low orbit): 278 m/s
    # Helioxis landing from orbit (same as low orbit to escape): 131 m/s
    required_dv = 28 + 271 + 278 + 131

    # Step 2: Generate all possible tank combinations
    tank_ids = list(tanks.keys())
    three_tank_combos = list(itertools.combinations(tank_ids, 3))
    four_tank_combos = list(itertools.combinations(tank_ids, 4))
    all_tank_combos = three_tank_combos + four_tank_combos

    min_cost = float('inf')
    best_config = None

    # Step 3 & 4: Iterate through all configurations and evaluate them
    for motor_id, motor_data in motors.items():
        for tank_combo in all_tank_combos:
            # Calculate total cost, wet mass (m0), and dry mass (mf)
            current_cost = motor_data['price']
            m0 = motor_data['mass']  # Start with motor mass
            mf = motor_data['mass']  # Start with motor mass

            for tank_id in tank_combo:
                current_cost += tanks[tank_id]['price']
                m0 += tanks[tank_id]['wet_mass']
                mf += tanks[tank_id]['dry_mass']

            # Calculate achievable Delta-V using Tsiolkovsky rocket equation
            # Δv = Ve * ln(m0 / mf)
            if mf > 0 and m0 > mf:
                achieved_dv = motor_data['ve'] * math.log(m0 / mf)
            else:
                achieved_dv = 0

            # Check if the rocket is viable and cheaper than the current best
            if achieved_dv >= required_dv:
                if current_cost < min_cost:
                    min_cost = current_cost
                    best_config = (motor_id, tank_combo)

    # Step 5: Format and print the final answer
    if best_config:
        motor_part = f"({best_config[0]})"
        tanks_part = ", ".join(map(str, sorted(best_config[1])))
        final_answer = f"{motor_part} {tanks_part}"
        print(final_answer)
        # The final answer is wrapped for the system to extract it.
        print(f"<<<{final_answer}>>>")
    else:
        print("No suitable rocket configuration found.")

solve_rocket_problem()