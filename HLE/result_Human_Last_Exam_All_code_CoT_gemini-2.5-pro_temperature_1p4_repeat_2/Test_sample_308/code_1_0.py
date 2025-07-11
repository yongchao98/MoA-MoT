import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """
    # --- Data from the tables ---
    motors = {
        1: {'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},  # ve converted to m/s
        2: {'price': 16000, 'mass': 1.54, 've': 2.3 * 1000} # ve converted to m/s
    }

    tanks = {
        1: {'price': 6000, 'wet': 5.2, 'dry': 3.9},
        2: {'price': 9000, 'wet': 7.8, 'dry': 5.1},
        3: {'price': 14000, 'wet': 11.1, 'dry': 6.0},
        4: {'price': 12000, 'wet': 10.1, 'dry': 7.5}
    }

    # --- Mission Delta-V Calculation ---
    # Delta-v Earth low orbit -> Helioxis low orbit
    dv_transfer = 271  # m/s
    # Delta-v Helioxis low orbit -> Helioxis surface
    dv_landing = 278 # m/s
    total_dv_required = dv_transfer + dv_landing

    valid_rockets = []
    tank_ids = list(tanks.keys())

    # --- Generate all tank combinations (3 or 4 tanks) ---
    tank_combos_3 = list(itertools.combinations(tank_ids, 3))
    tank_combos_4 = list(itertools.combinations(tank_ids, 4))
    all_tank_combos = tank_combos_3 + tank_combos_4

    # --- Iterate through all possible rocket configurations ---
    for motor_id, motor_data in motors.items():
        for combo in all_tank_combos:
            # Calculate totals for this configuration
            current_price = motor_data['price']
            total_wet_mass = motor_data['mass']
            total_dry_mass = motor_data['mass']

            for tank_id in combo:
                tank_data = tanks[tank_id]
                current_price += tank_data['price']
                total_wet_mass += tank_data['wet']
                total_dry_mass += tank_data['dry']

            # Calculate rocket's delta-v using Tsiolkovsky rocket equation
            mass_ratio = total_wet_mass / total_dry_mass
            rocket_dv = motor_data['ve'] * math.log(mass_ratio)

            # Check if the rocket meets the delta-v requirement
            if rocket_dv >= total_dv_required:
                valid_rockets.append({
                    'motor': motor_id,
                    'tanks': sorted(combo),
                    'price': current_price,
                    'dv': rocket_dv
                })

    # --- Find the cheapest rocket among the valid ones ---
    if not valid_rockets:
        print("No rocket configuration meets the mission's delta-v requirement.")
    else:
        cheapest_rocket = min(valid_rockets, key=lambda x: x['price'])
        
        motor_num = cheapest_rocket['motor']
        tank_nums_str = ", ".join(map(str, cheapest_rocket['tanks']))
        
        # Format and print the final answer
        result_string = f"({motor_num}) {tank_nums_str}"
        print(result_string)
        
        # Also enclose in the final answer format for the platform
        print(f"<<<{result_string}>>>")

if __name__ == '__main__':
    find_cheapest_rocket()