import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a given mission based on provided component data.
    """
    # --- Data Definition ---
    # Motor data: {id: {'price': price, 'mass': mass, 've': exhaust_velocity_m/s}}
    motors = {
        1: {'price': 8000, 'mass': 1.3, 've': 1200},
        2: {'price': 16000, 'mass': 1.54, 've': 2300}
    }

    # Tank data: {id: {'price': price, 'wet_mass': wet_mass, 'dry_mass': dry_mass}}
    tanks = {
        1: {'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        2: {'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        3: {'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        4: {'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    }

    # --- Mission Requirements ---
    delta_v_transfer = 271  # Earth to Helioxis
    delta_v_landing = 278   # Helioxis orbit to surface
    required_delta_v = delta_v_transfer + delta_v_landing

    valid_configs = []

    # --- Iterate through all possible configurations ---
    # Iterate through each motor
    for motor_id, motor_data in motors.items():
        # Iterate through allowed number of tanks (3 or 4)
        for num_tanks in [3, 4]:
            # Generate all combinations of tanks of the specified size
            for tank_combo_ids in itertools.combinations(tanks.keys(), num_tanks):
                
                # Calculate properties for the current configuration
                current_cost = motor_data['price']
                m_initial = motor_data['mass']  # Start with motor mass
                m_final = motor_data['mass']    # Motor mass is the same for wet and dry

                for tank_id in tank_combo_ids:
                    current_cost += tanks[tank_id]['price']
                    m_initial += tanks[tank_id]['wet_mass']
                    m_final += tanks[tank_id]['dry_mass']

                # Calculate achievable delta-v using Tsiolkovsky rocket equation
                # dV = Ve * ln(M_initial / M_final)
                achieved_delta_v = motor_data['ve'] * math.log(m_initial / m_final)

                # Check if the configuration meets the mission requirement
                if achieved_delta_v >= required_delta_v:
                    valid_configs.append({
                        'motor': motor_id,
                        'tanks': sorted(tank_combo_ids),
                        'cost': current_cost,
                        'delta_v': achieved_delta_v
                    })

    # --- Find the cheapest valid configuration ---
    if not valid_configs:
        print("No valid configuration found that meets the delta-v requirement.")
        return

    cheapest_config = min(valid_configs, key=lambda x: x['cost'])

    # --- Output the results ---
    print("The required delta-v for the mission is calculated as:")
    print(f"{delta_v_transfer} + {delta_v_landing} = {required_delta_v}")
    
    motor_part = f"({cheapest_config['motor']})"
    tanks_part = ", ".join(map(str, cheapest_config['tanks']))
    final_answer_string = f"{motor_part} {tanks_part}"
    
    print("\nThe cheapest valid configuration is:")
    print(final_answer_string)
    
    print(f"\n<<<{final_answer_string}>>>")


find_cheapest_rocket()