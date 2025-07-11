import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a trip from low Earth orbit
    to the surface of Helioxis based on the provided component and travel data.
    """
    # Step 1: Define component data and calculate required Delta-v
    motors = [
        {'id': 1, 'price': 8000, 'mass': 1.3, 've': 1200},  # ve is pre-converted to m/s
        {'id': 2, 'price': 16000, 'mass': 1.54, 've': 2300}
    ]

    tanks = [
        {'id': 1, 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # Delta-v from low Earth orbit to low Helioxis orbit
    dv_transfer = 271  # m/s
    # Delta-v from low Helioxis orbit to surface (same as surface to low orbit)
    dv_landing = 278  # m/s
    total_dv_required = dv_transfer + dv_landing

    # Step 2: Enumerate all configurations
    valid_configs = []
    tank_indices = list(range(len(tanks)))
    
    # Generate combinations for 3 and 4 tanks
    combos_3_tanks = itertools.combinations(tank_indices, 3)
    combos_4_tanks = itertools.combinations(tank_indices, 4)
    all_tank_combos = list(combos_3_tanks) + list(combos_4_tanks)

    # Step 3 & 4: Calculate performance/cost for each config and filter
    for motor in motors:
        for tank_combo_indices in all_tank_combos:
            current_tanks = [tanks[i] for i in tank_combo_indices]

            # Calculate cost
            config_cost = motor['price'] + sum(t['price'] for t in current_tanks)

            # Calculate masses
            total_wet_mass = motor['mass'] + sum(t['wet_mass'] for t in current_tanks)
            total_dry_mass = motor['mass'] + sum(t['dry_mass'] for t in current_tanks)

            # Calculate achievable delta-v using Tsiolkovsky rocket equation
            if total_dry_mass > 0 and total_wet_mass > total_dry_mass:
                mass_ratio = total_wet_mass / total_dry_mass
                achievable_dv = motor['ve'] * math.log(mass_ratio)
            else:
                achievable_dv = 0

            # Check if configuration meets the delta-v requirement
            if achievable_dv >= total_dv_required:
                tank_ids = sorted([t['id'] for t in current_tanks])
                valid_configs.append({
                    'motor_id': motor['id'],
                    'tank_ids': tank_ids,
                    'cost': config_cost
                })

    # Step 5 & 6: Select the cheapest option and format the output
    if not valid_configs:
        print("No valid rocket configuration found to complete the mission.")
    else:
        cheapest_config = min(valid_configs, key=lambda x: x['cost'])
        
        motor_id_str = str(cheapest_config['motor_id'])
        tank_ids_str = ", ".join(map(str, cheapest_config['tank_ids']))
        
        # The final answer format is "(motor) tank1, tank2, tank3"
        final_answer = f"({motor_id_str}) {tank_ids_str}"
        print(final_answer)

find_cheapest_rocket()