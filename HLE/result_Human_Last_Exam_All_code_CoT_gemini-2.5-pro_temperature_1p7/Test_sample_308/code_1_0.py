import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """
    # --- Data from the tables ---
    motors = [
        {'id': 1, 'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},
        {'id': 2, 'price': 16000, 'mass': 1.54, 've': 2.3 * 1000}
    ]

    tanks = [
        {'id': 1, 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- 1. Calculate Required Delta-V ---
    # Delta-v values in m/s
    delta_v_earth_escape = 28       # Earth: low orbit to escape
    delta_v_transfer = 271          # Earth to Helioxis
    delta_v_helioxis_landing = 278  # Helioxis: surface to low orbit

    total_required_delta_v = delta_v_earth_escape + delta_v_transfer + delta_v_helioxis_landing

    # --- Initialize variables to track the best configuration ---
    min_cost = float('inf')
    best_config = None

    # --- 2. Generate all tank combinations (3 or 4 tanks) ---
    tank_combinations_3 = list(itertools.combinations(tanks, 3))
    tank_combinations_4 = list(itertools.combinations(tanks, 4))
    all_tank_combos = tank_combinations_3 + tank_combinations_4
    
    # --- 3. Iterate through all configurations and calculate performance ---
    for motor in motors:
        for tank_combo in all_tank_combos:
            # Calculate cost for the current configuration
            current_cost = motor['price'] + sum(t['price'] for t in tank_combo)

            # Calculate mass for the Tsiolkovsky equation
            # M0 = initial (wet) mass, Mf = final (dry) mass
            M0 = motor['mass'] + sum(t['wet_mass'] for t in tank_combo)
            Mf = motor['mass'] + sum(t['dry_mass'] for t in tank_combo)

            # Get exhaust velocity (Ve)
            Ve = motor['ve']

            # Calculate achievable delta-v
            achieved_delta_v = Ve * math.log(M0 / Mf)

            # --- 4. Identify the optimal rocket ---
            # Check if it meets the delta-v requirement and is cheaper than the best one found so far
            if achieved_delta_v >= total_required_delta_v and current_cost < min_cost:
                min_cost = current_cost
                best_config = {
                    'motor': motor,
                    'tanks': tank_combo,
                    'cost': current_cost,
                    'delta_v': achieved_delta_v
                }

    # --- 5. Format and print the output ---
    if best_config:
        motor_id = best_config['motor']['id']
        # Sort tank IDs for consistent output
        tank_ids = sorted([t['id'] for t in best_config['tanks']])
        tank_ids_str = ', '.join(map(str, tank_ids))
        
        final_answer_string = f"({motor_id}) {tank_ids_str}"
        
        print(f"The cheapest rocket capable of traveling from low Earth orbit to the surface of Helioxis is: {final_answer_string}")
        print(f"<<<{final_answer_string}>>>")
    else:
        print("No suitable rocket configuration was found that meets the mission requirements.")

if __name__ == '__main__':
    find_cheapest_rocket()