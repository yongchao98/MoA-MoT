import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis.
    """
    # --- Data from the tables ---
    motors = [
        {'id': 1, 'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 'exhaust_velocity': 1200},
        {'id': 2, 'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 'exhaust_velocity': 2300}
    ]

    tanks = [
        {'id': 1, 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- 1. Calculate the total Delta-v requirement ---
    dv_earth_escape = 28  # Delta-v (low orbit to escape) for Earth
    dv_transfer = 271      # Optimal transfer delta-v from Earth to Helioxis
    dv_helioxis_landing = 278 # Delta-v (surface to low orbit) for Helioxis
    
    total_dv_required = dv_earth_escape + dv_transfer + dv_helioxis_landing
    
    print("Required Delta-v Calculation:")
    print(f"{dv_earth_escape} (Earth escape) + {dv_transfer} (Transfer) + {dv_helioxis_landing} (Helioxis landing) = {total_dv_required} m/s")
    print("-" * 30)

    valid_configurations = []

    # --- 2. Generate and test all configurations ---
    for motor in motors:
        # Combinations of 3 or 4 tanks
        for num_tanks in [3, 4]:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                # 3. Calculate performance
                current_cost = motor['price'] + sum(t['price'] for t in tank_combo)
                
                m_initial = motor['wet_mass'] + sum(t['wet_mass'] for t in tank_combo)
                m_final = motor['dry_mass'] + sum(t['dry_mass'] for t in tank_combo)
                
                v_e = motor['exhaust_velocity']

                # Tsiolkovsky rocket equation
                if m_initial > m_final:
                    delta_v_achieved = v_e * math.log(m_initial / m_final)
                else:
                    delta_v_achieved = 0

                # 4. Identify valid configurations
                if delta_v_achieved >= total_dv_required:
                    tank_ids = sorted([t['id'] for t in tank_combo])
                    valid_configurations.append({
                        'cost': current_cost,
                        'motor_id': motor['id'],
                        'tank_ids': tank_ids
                    })

    # --- 5. Find the cheapest valid configuration and print ---
    if not valid_configurations:
        print("No suitable rocket configuration found.")
    else:
        cheapest_config = min(valid_configurations, key=lambda x: x['cost'])
        
        motor_id = cheapest_config['motor_id']
        tank_ids_str = ", ".join(map(str, cheapest_config['tank_ids']))
        
        print("Cheapest valid rocket configuration:")
        print(f"({motor_id}) {tank_ids_str}")

if __name__ == '__main__':
    find_cheapest_rocket()