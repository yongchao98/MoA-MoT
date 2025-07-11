import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a given mission based on provided component data.
    """
    # --- Step 1: Define Data and Mission Requirements ---

    # Component data from tables
    motors = {
        1: {'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 've': 1200},  # ve converted to m/s
        2: {'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 've': 2300}
    }

    tanks = {
        1: {'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        2: {'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        3: {'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        4: {'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    }

    # Mission: Low Earth orbit to the surface of Helioxis
    # From tables:
    # - Optimal transfer delta-v (Earth to Helioxis) = 271 m/s
    # - Delta-v (surface to low orbit, Helioxis) = 278 m/s (landing is the reverse)
    required_delta_v = 271 + 278

    # --- Step 2 & 3: Iterate through configurations and calculate performance ---

    min_cost = float('inf')
    best_config = None

    # Iterate through each motor
    for motor_id, motor_data in motors.items():
        # Iterate through combination sizes (3 or 4 tanks)
        for num_tanks in [3, 4]:
            # Generate all combinations of tanks for the given size
            for tank_combo_ids in itertools.combinations(tanks.keys(), num_tanks):
                
                # Calculate total properties for the current configuration
                current_cost = motor_data['price']
                total_wet_mass = motor_data['wet_mass']
                total_dry_mass = motor_data['dry_mass']
                
                for tank_id in tank_combo_ids:
                    current_cost += tanks[tank_id]['price']
                    total_wet_mass += tanks[tank_id]['wet_mass']
                    total_dry_mass += tanks[tank_id]['dry_mass']
                
                # Calculate delta-v using Tsiolkovsky rocket equation
                # Δv = Ve * ln(M_initial / M_final)
                if total_dry_mass > 0:
                    mass_ratio = total_wet_mass / total_dry_mass
                    # The rocket only has delta-v if it has fuel (mass_ratio > 1)
                    delta_v = motor_data['ve'] * math.log(mass_ratio) if mass_ratio > 1 else 0
                else:
                    delta_v = 0

                # --- Step 4: Check if the configuration is successful and cheapest ---
                if delta_v >= required_delta_v:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        best_config = {
                            'motor': motor_id,
                            'tanks': sorted(tank_combo_ids)
                        }

    # --- Step 5: Format and print the final answer ---
    if best_config:
        motor_num = best_config['motor']
        tank_nums_str = ", ".join(map(str, best_config['tanks']))
        result_string = f"({motor_num}) {tank_nums_str}"
        print(result_string)
    else:
        print("No suitable rocket configuration was found.")

solve_rocket_problem()
<<<(2) 1, 2, 4>>>