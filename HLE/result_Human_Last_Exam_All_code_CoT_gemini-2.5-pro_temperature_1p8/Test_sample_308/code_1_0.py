import math
import itertools

def solve():
    """
    Finds the cheapest rocket configuration for a trip from low Earth orbit to the surface of Helioxis.
    """
    # Step 1: Define required delta-v for the mission
    # Delta-v (Earth orbit to Helioxis orbit) + Delta-v (Helioxis orbit to surface)
    required_dv = 271 + 278  # in m/s

    # Define component data
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'mass': 1.3, 'v_e': 1.2 * 1000},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'mass': 1.54, 'v_e': 2.3 * 1000}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    cheapest_config = None
    min_cost = float('inf')

    # Step 2: Iterate through all possible valid configurations
    for motor in motors:
        # Combinations of 3 or 4 tanks
        for num_tanks in [3, 4]:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                
                # Calculate total cost
                current_cost = motor['price'] + sum(t['price'] for t in tank_combo)

                # Optimization: if it's already more expensive than the best we've found, skip it
                if current_cost >= min_cost:
                    continue

                # Calculate total wet mass (m0) and dry mass (mf)
                m0 = motor['mass'] + sum(t['wet_mass'] for t in tank_combo)
                mf = motor['mass'] + sum(t['dry_mass'] for t in tank_combo)

                # Step 3: Calculate the rocket's delta-v using Tsiolkovsky equation
                # Handle cases where mf might be zero or m0/mf <= 1 to avoid math domain errors
                if mf > 0 and m0 > mf:
                    achieved_dv = motor['v_e'] * math.log(m0 / mf)
                else:
                    achieved_dv = 0

                # Step 4 & 5: Check if the rocket is capable and the cheapest so far
                if achieved_dv >= required_dv:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        tank_ids = sorted([t['id'] for t in tank_combo])
                        cheapest_config = {
                            'motor_id': motor['id'],
                            'tank_ids': tank_ids,
                            'cost': current_cost,
                            'dv': achieved_dv
                        }

    # Step 6: Format and print the output
    if cheapest_config:
        motor_id = cheapest_config['motor_id']
        tank_ids_str = ", ".join(map(str, cheapest_config['tank_ids']))
        result_str = f"({motor_id}) {tank_ids_str}"
        print(result_str)
        # Final answer block
        print(f"<<<{result_str}>>>")
    else:
        print("No capable rocket configuration found.")

solve()