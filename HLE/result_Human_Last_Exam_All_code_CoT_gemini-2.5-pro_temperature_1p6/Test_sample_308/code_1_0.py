import math
import itertools

def solve_rocket_problem():
    """
    Finds the cheapest rocket configuration (1 motor, 3 or 4 tanks)
    to travel from low Earth orbit to the surface of Helioxis.
    """
    # --- Data Setup ---
    # Store motor and tank data. Using dictionaries for clarity.
    # Masses are in metric tons (mt), exhaust velocity in km/s.
    motors = [
        {'id': 1, 'price': 8000, 'mass': 1.3, 've_kms': 1.2},
        {'id': 2, 'price': 16000, 'mass': 1.54, 've_kms': 2.3}
    ]
    tanks = [
        {'id': 1, 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- Step 1: Calculate Delta-V Requirement ---
    # Δv from Earth orbit to Helioxis orbit + Δv from Helioxis orbit to surface
    delta_v_required_ms = 271 + 278

    # --- Initialize variables to track the best solution ---
    cheapest_capable_config = None
    min_cost = float('inf')

    # --- Step 2 & 3: Iterate through all configurations and calculate performance ---
    for motor in motors:
        # The constraints allow for combinations of 3 or 4 tanks
        for num_tanks in [3, 4]:
            # Generate all combinations of tanks for the current size
            for tank_combo in itertools.combinations(tanks, num_tanks):
                
                # a. Calculate total cost for the current configuration
                current_cost = motor['price'] + sum(t['price'] for t in tank_combo)

                # b. Calculate initial and final masses (in metric tons)
                initial_mass = motor['mass'] + sum(t['wet_mass'] for t in tank_combo)
                final_mass = motor['mass'] + sum(t['dry_mass'] for t in tank_combo)

                # c. Calculate achievable delta-v using Tsiolkovsky rocket equation
                # Convert exhaust velocity from km/s to m/s
                ve_ms = motor['ve_kms'] * 1000
                
                # Mass ratio must be > 1 for positive delta-v
                if initial_mass > final_mass:
                    mass_ratio = initial_mass / final_mass
                    achieved_delta_v = ve_ms * math.log(mass_ratio)
                else:
                    achieved_delta_v = 0
                
                # --- Step 4: Check if configuration is capable and the cheapest ---
                if achieved_delta_v >= delta_v_required_ms:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        # Store the details of this new best configuration
                        tank_ids = sorted([t['id'] for t in tank_combo])
                        cheapest_capable_config = {
                            'motor_id': motor['id'],
                            'tank_ids': tank_ids,
                            'cost': current_cost
                        }

    # --- Output the Final Answer ---
    if cheapest_capable_config:
        motor_id = cheapest_capable_config['motor_id']
        tank_ids_str = ", ".join(map(str, cheapest_capable_config['tank_ids']))
        # Print the final result in the format: (motor_id) tank1_id, tank2_id, ...
        print(f"({motor_id}) {tank_ids_str}")
    else:
        print("No capable rocket configuration found within the given constraints.")

solve_rocket_problem()