import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a given space mission.
    """
    # Step 1: Define the component data and mission requirements.
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'mass': 1.3, 've': 1200},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'mass': 1.54, 've': 2300}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # Delta-v requirement for Earth orbit -> Helioxis orbit -> Helioxis surface
    delta_v_transfer = 271  # m/s
    delta_v_landing = 278   # m/s
    delta_v_required = delta_v_transfer + delta_v_landing

    valid_configs = []

    # Step 2 & 3: Generate all configurations and calculate their performance.
    # Iterate through motors
    for motor in motors:
        # Iterate through combinations of 3 and 4 tanks
        for num_tanks in [3, 4]:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                
                # Calculate total cost
                config_price = motor['price'] + sum(t['price'] for t in tank_combo)
                
                # Calculate initial (m0) and final (mf) mass
                m0 = motor['mass'] + sum(t['wet_mass'] for t in tank_combo)
                mf = motor['mass'] + sum(t['dry_mass'] for t in tank_combo)
                
                # Calculate achievable delta-v using the Tsiolkovsky rocket equation
                if m0 > mf:
                    achieved_delta_v = motor['ve'] * math.log(m0 / mf)
                else:
                    achieved_delta_v = 0
                
                # Step 4: Check if the configuration is valid
                if achieved_delta_v >= delta_v_required:
                    tank_ids = sorted([t['id'] for t in tank_combo])
                    valid_configs.append({
                        'price': config_price,
                        'motor_id': motor['id'],
                        'tank_ids': tank_ids
                    })

    # Step 5: Find the cheapest valid configuration.
    if not valid_configs:
        print(f"The total delta-v requirement is {delta_v_required} m/s.")
        print("No valid configuration found to meet the delta-v requirement.")
    else:
        cheapest_config = min(valid_configs, key=lambda x: x['price'])
        
        # Step 6: Format and print the final output.
        motor_id_str = f"({cheapest_config['motor_id']})"
        tank_ids_str = ', '.join(map(str, cheapest_config['tank_ids']))
        final_answer_str = f"{motor_id_str} {tank_ids_str};"
        
        print(f"The total delta-v requirement is the sum of delta-v for orbital transfer and landing.")
        print(f"{delta_v_transfer} + {delta_v_landing} = {delta_v_required}")
        print(f"The cheapest rocket configuration capable of the mission is: {final_answer_str}")
        print(f"<<<{final_answer_str}>>>")

if __name__ == '__main__':
    find_cheapest_rocket()