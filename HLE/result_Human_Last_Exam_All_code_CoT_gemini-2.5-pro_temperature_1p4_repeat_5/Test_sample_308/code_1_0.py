import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a trip from low Earth orbit
    to the surface of Helioxis based on the provided component data and constraints.
    """

    # --- Step 1: Define data and calculate required Delta-v ---

    # Component Data
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

    # Delta-v requirements (in m/s)
    dv_leo_to_escape = 28        # Earth: low orbit to escape
    dv_earth_to_helioxis = 271   # Optimal transfer Earth -> Helioxis
    dv_helioxis_orbit_to_surface = 278 # Helioxis: surface to low orbit

    required_dv = dv_leo_to_escape + dv_earth_to_helioxis + dv_helioxis_orbit_to_surface

    # --- Step 2 & 3: Iterate through configurations and evaluate them ---

    capable_rockets = []
    
    # Iterate through each motor
    for motor in motors:
        # Iterate through combinations of 3 and 4 tanks
        for num_tanks in [3, 4]:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                # Calculate total properties for this configuration
                total_price = motor['price'] + sum(t['price'] for t in tank_combo)
                total_wet_mass = motor['wet_mass'] + sum(t['wet_mass'] for t in tank_combo)
                total_dry_mass = motor['dry_mass'] + sum(t['dry_mass'] for t in tank_combo)
                
                # Tsiolkovsky rocket equation
                if total_dry_mass > 0 and total_wet_mass > total_dry_mass:
                    achieved_dv = motor['ve'] * math.log(total_wet_mass / total_dry_mass)
                else:
                    achieved_dv = 0

                # Check if this configuration is capable
                if achieved_dv >= required_dv:
                    tank_ids = sorted([t['id'] for t in tank_combo])
                    capable_rockets.append({
                        'motor_id': motor['id'],
                        'tank_ids': tank_ids,
                        'price': total_price,
                        'delta_v': achieved_dv
                    })

    # --- Step 4: Find the cheapest solution among capable rockets ---

    if not capable_rockets:
        print("No capable rocket configuration found.")
        return

    cheapest_rocket = min(capable_rockets, key=lambda x: x['price'])
    
    # --- Step 5: Format and print the final output ---
    
    motor_id_str = str(cheapest_rocket['motor_id'])
    tank_ids_str = ', '.join(map(str, cheapest_rocket['tank_ids']))
    
    final_answer = f"({motor_id_str}) {tank_ids_str}"
    
    print(final_answer)
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    find_cheapest_rocket()