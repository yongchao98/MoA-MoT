import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a given space mission.
    """

    # Data from the tables
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 've': 1.2 * 1000},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 've': 2.3 * 1000}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # Required delta-v calculation
    # Earth low orbit to Helioxis low orbit = 271 m/s
    # Helioxis low orbit to surface = 278 m/s
    required_dv = 271 + 278

    # Generate all possible tank combinations (3 or 4 tanks)
    tank_combos_3 = list(itertools.combinations(tanks, 3))
    tank_combos_4 = list(itertools.combinations(tanks, 4))
    all_tank_combos = tank_combos_3 + tank_combos_4

    capable_rockets = []

    # Iterate through all possible configurations
    for motor in motors:
        for tank_combo in all_tank_combos:
            # Calculate total cost, wet mass, and dry mass for the configuration
            config_cost = motor['price'] + sum(tank['price'] for tank in tank_combo)
            m_initial = motor['wet_mass'] + sum(tank['wet_mass'] for tank in tank_combo)
            m_final = motor['dry_mass'] + sum(tank['dry_mass'] for tank in tank_combo)
            
            # Tsiolkovsky rocket equation
            if m_initial > m_final:
                achieved_dv = motor['ve'] * math.log(m_initial / m_final)
            else:
                achieved_dv = 0

            # Check if the rocket is capable
            if achieved_dv >= required_dv:
                capable_rockets.append({
                    'motor_id': motor['id'],
                    'tank_ids': [tank['id'] for tank in tank_combo],
                    'cost': config_cost,
                    'achieved_dv': achieved_dv,
                    've': motor['ve'],
                    'm_initial': m_initial,
                    'm_final': m_final
                })

    # Find the cheapest capable rocket
    if not capable_rockets:
        print("No capable rocket configuration found.")
        return

    cheapest_rocket = min(capable_rockets, key=lambda x: x['cost'])

    # Format and print the final answer
    motor_num = cheapest_rocket['motor_id']
    tank_nums = sorted(cheapest_rocket['tank_ids'])
    tank_str = ", ".join(map(str, tank_nums))
    answer_str = f"({motor_num}) {tank_str};"
    
    ve = cheapest_rocket['ve']
    m_initial = round(cheapest_rocket['m_initial'], 2)
    m_final = round(cheapest_rocket['m_final'], 2)
    dv = round(cheapest_rocket['achieved_dv'], 1)

    equation_str = f"{dv} = {int(ve)} * ln({m_initial} / {m_final})"

    print(answer_str)
    print("The final equation is:")
    print(equation_str)
    print(f"<<<({motor_num}) {tank_str};>>>")


solve_rocket_problem()