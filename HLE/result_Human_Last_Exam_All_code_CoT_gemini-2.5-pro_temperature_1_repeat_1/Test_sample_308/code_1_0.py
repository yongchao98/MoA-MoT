import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a journey from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """

    # --- Data Setup ---
    # Store component data in dictionaries for easy access.
    # Motor IDs are 1 and 2. Tank IDs are 1, 2, 3, and 4.
    motors = {
        1: {'price': 8000, 'mass': 1.3, 'v_e': 1.2 * 1000},  # Convert v_e to m/s
        2: {'price': 16000, 'mass': 1.54, 'v_e': 2.3 * 1000} # Convert v_e to m/s
    }

    tanks = {
        1: {'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        2: {'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        3: {'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        4: {'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    }
    tank_ids = list(tanks.keys())

    # --- Mission Requirements ---
    # Delta-v from Earth (low orbit) to Helioxis (low orbit transfer)
    dv_transfer = 271  # m/s
    # Delta-v from Helioxis (low orbit) to Helioxis (surface landing)
    dv_landing = 278  # m/s
    # Total delta-v required for the mission
    required_dv = dv_transfer + dv_landing

    # --- Analysis ---
    min_cost = float('inf')
    best_config = None

    # Iterate through each motor option
    for motor_id, motor_data in motors.items():
        
        # Consider combinations of 3 and 4 tanks, as per the constraint
        for num_tanks in [3, 4]:
            
            # Generate all unique combinations of tanks for the given number
            for tank_combo_ids in itertools.combinations(tank_ids, num_tanks):
                
                # --- Calculate properties for the current rocket configuration ---
                
                # 1. Total cost
                current_cost = motor_data['price'] + sum(tanks[tid]['price'] for tid in tank_combo_ids)
                
                # 2. Total initial (wet) mass (m0)
                m0 = motor_data['mass'] + sum(tanks[tid]['wet_mass'] for tid in tank_combo_ids)
            
                # 3. Total final (dry) mass (mf)
                mf = motor_data['mass'] + sum(tanks[tid]['dry_mass'] for tid in tank_combo_ids)
            
                # 4. Exhaust velocity (ve)
                ve = motor_data['v_e']
                
                # --- Tsiolkovsky Rocket Equation: delta_v = ve * ln(m0 / mf) ---
                # Check for valid mass ratio to avoid math errors
                if m0 > mf > 0:
                    achieved_dv = ve * math.log(m0 / mf)
                else:
                    achieved_dv = 0

                # --- Check if this configuration is a new best solution ---
                # It must meet the delta-v requirement and be cheaper than any previous solution
                if achieved_dv >= required_dv:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        best_config = {
                            'motor': motor_id,
                            'tanks': sorted(list(tank_combo_ids)),
                        }

    # --- Output the Final Answer ---
    if best_config:
        # Format the tank list into a "1, 2, 3" string
        tanks_str = ", ".join(map(str, best_config['tanks']))
        # Construct the final answer string in the required format
        final_answer = f"({best_config['motor']}) {tanks_str}"
        print(final_answer)
    else:
        # This will be printed if no combination meets the delta-v requirement
        print("No suitable configuration found.")

solve_rocket_problem()