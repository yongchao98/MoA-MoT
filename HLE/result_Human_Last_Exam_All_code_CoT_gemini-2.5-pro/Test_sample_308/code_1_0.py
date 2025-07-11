import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a given space mission.
    """

    # --- Data from tables ---
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'mass': 1.3, 've_kms': 1.2},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'mass': 1.54, 've_kms': 2.3}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- Mission Requirements ---
    # Delta-v from low Earth orbit to low Helioxis orbit
    dv_transfer = 271  # m/s
    # Delta-v from low Helioxis orbit to surface
    dv_landing = 278  # m/s
    required_delta_v = dv_transfer + dv_landing

    # --- Generate all possible configurations ---
    tank_ids = [tank['id'] for tank in tanks]
    three_tank_combos = list(itertools.combinations(tank_ids, 3))
    four_tank_combos = list(itertools.combinations(tank_ids, 4))
    all_tank_combos = three_tank_combos + four_tank_combos

    valid_configurations = []

    # --- Evaluate each configuration ---
    for motor in motors:
        for tank_combo_ids in all_tank_combos:
            selected_tanks = [tank for tank in tanks if tank['id'] in tank_combo_ids]

            # Calculate total cost, masses
            total_price = motor['price'] + sum(tank['price'] for tank in selected_tanks)
            total_wet_mass_tanks = sum(tank['wet_mass'] for tank in selected_tanks)
            total_dry_mass_tanks = sum(tank['dry_mass'] for tank in selected_tanks)

            m_initial = motor['mass'] + total_wet_mass_tanks
            m_final = motor['mass'] + total_dry_mass_tanks
            
            # Exhaust velocity in m/s
            ve_ms = motor['ve_kms'] * 1000

            # Calculate delta-v using Tsiolkovsky Rocket Equation
            # Check for m_final > 0 to avoid division by zero
            if m_final > 0 and m_initial > m_final:
                delta_v = ve_ms * math.log(m_initial / m_final)
            else:
                delta_v = 0

            # Check if configuration meets the delta-v requirement
            if delta_v >= required_delta_v:
                valid_configurations.append({
                    'motor_id': motor['id'],
                    'tank_ids': tank_combo_ids,
                    'price': total_price,
                    'delta_v': delta_v
                })

    # --- Find the cheapest valid configuration ---
    if not valid_configurations:
        print("No valid rocket configuration found.")
        return

    cheapest_config = min(valid_configurations, key=lambda x: x['price'])
    
    motor_id = cheapest_config['motor_id']
    tank_ids_str = ', '.join(map(str, sorted(cheapest_config['tank_ids'])))
    
    # Print the final answer in the required format
    print(f"({motor_id}) {tank_ids_str}")

# Run the solver
solve_rocket_problem()