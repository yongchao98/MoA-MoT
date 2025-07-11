import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis, subject to component and assembly constraints.
    """
    # Define rocket component data
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

    # 1. Determine Mission Delta-V
    # Optimal transfer delta-v from Earth to Helioxis
    dv_transfer = 271  # m/s
    # Delta-v from Helioxis surface to low orbit (same as low orbit to surface)
    dv_landing = 278  # m/s
    required_dv = dv_transfer + dv_landing

    # 3. Iterate Through All Configurations
    tank_ids = list(tanks.keys())
    combinations_3_tanks = list(itertools.combinations(tank_ids, 3))
    combinations_4_tanks = list(itertools.combinations(tank_ids, 4))
    all_tank_combinations = combinations_3_tanks + combinations_4_tanks

    valid_configs = []

    for motor_id, motor_data in motors.items():
        for tank_combo in all_tank_combinations:
            # 4. Evaluate Each Configuration
            config_cost = motor_data['price']
            m_initial = motor_data['mass']
            m_final = motor_data['mass']

            for tank_id in tank_combo:
                config_cost += tanks[tank_id]['price']
                m_initial += tanks[tank_id]['wet_mass']
                m_final += tanks[tank_id]['dry_mass']

            ve = motor_data['ve']

            # 2. Implement the Rocket Equation
            achieved_dv = ve * math.log(m_initial / m_final)

            # Check if the configuration is viable
            if achieved_dv >= required_dv:
                valid_configs.append({
                    'cost': config_cost,
                    'motor': motor_id,
                    'tanks': sorted(tank_combo),
                    'dv': achieved_dv,
                    've': ve,
                    'm_initial': m_initial,
                    'm_final': m_final
                })

    # 5. Find the Cheapest Viable Rocket
    if not valid_configs:
        print("No rocket configuration meets the mission requirements.")
        return

    cheapest_config = min(valid_configs, key=lambda x: x['cost'])

    # 6. Output the Result
    motor = cheapest_config['motor']
    tanks_list = cheapest_config['tanks']
    cost = cheapest_config['cost']
    dv = cheapest_config['dv']
    ve = cheapest_config['ve']
    m_initial = cheapest_config['m_initial']
    m_final = cheapest_config['m_final']
    
    tanks_str = ", ".join(map(str, tanks_list))

    print(f"The cheapest rocket configuration is Motor {motor} with tanks {tanks_str}.")
    print(f"The required delta-v is {required_dv} m/s, and this rocket achieves {dv:.1f} m/s.")
    print("\nThe Tsiolkovsky rocket equation for this configuration, with all numbers, is:")
    print(f"Delta-v = {ve} * ln({m_initial:.2f} / {m_final:.2f})")
    
    final_answer_format = f"({motor}) {tanks_str}"
    print(f"\nFinal Answer: {final_answer_format}")

if __name__ == '__main__':
    find_cheapest_rocket()