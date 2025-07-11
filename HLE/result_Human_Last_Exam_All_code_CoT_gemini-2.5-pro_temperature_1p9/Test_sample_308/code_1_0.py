import math
import itertools

def find_cheapest_rocket():
    """
    This script determines the cheapest rocket configuration (one motor and 3 or 4 tanks)
    capable of traveling from low Earth orbit to the surface of Helioxis.
    """
    # --- Data from tables ---
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 'exhaust_velocity': 1.2},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 'exhaust_velocity': 2.3}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- Step 1: Calculate Required Delta-V ---
    # Transfer from low Earth orbit to low Helioxis orbit
    dv_transfer_earth_to_helioxis = 271  # m/s
    # Landing on Helioxis from low orbit (same as surface to low orbit)
    dv_landing_on_helioxis = 278 # m/s
    total_required_dv = dv_transfer_earth_to_helioxis + dv_landing_on_helioxis

    # --- Step 2: Iterate through configurations to find the cheapest valid one ---
    min_cost = float('inf')
    best_config = None

    # Generate combinations of 3 and 4 tanks
    tank_indices = list(range(len(tanks)))
    tank_combos_3 = itertools.combinations(tank_indices, 3)
    tank_combos_4 = itertools.combinations(tank_indices, 4)
    all_tank_combos = list(tank_combos_3) + list(tank_combos_4)

    for motor in motors:
        for tank_index_combo in all_tank_combos:
            current_tanks = [tanks[i] for i in tank_index_combo]
            
            # --- Calculate properties for the current configuration ---
            total_cost = motor['price'] + sum(t['price'] for t in current_tanks)
            
            # Tsiolkovsky Rocket Equation inputs
            m0 = motor['wet_mass'] + sum(t['wet_mass'] for t in current_tanks) # Total wet mass
            mf = motor['dry_mass'] + sum(t['dry_mass'] for t in current_tanks) # Total dry mass
            ve = motor['exhaust_velocity'] * 1000 # Exhaust velocity in m/s
            
            # Ensure no math errors if mf is greater or equal to m0
            if m0 > mf:
                calculated_dv = ve * math.log(m0 / mf)
            else:
                calculated_dv = 0

            # --- Check if this configuration is a valid and cheaper candidate ---
            if calculated_dv >= total_required_dv:
                if total_cost < min_cost:
                    min_cost = total_cost
                    best_config = {
                        'motor': motor,
                        'tanks': current_tanks,
                        'cost': total_cost,
                        'calculated_dv': calculated_dv,
                        've': ve,
                        'm0': m0,
                        'mf': mf
                    }

    # --- Step 3: Print the final results ---
    if best_config:
        motor_id = best_config['motor']['id']
        tank_ids = sorted([t['id'] for t in best_config['tanks']])
        tank_ids_str = ', '.join(map(str, tank_ids))
        
        print(f"The cheapest rocket capable of the mission is composed of motor {motor_id} and tanks {tank_ids_str}.")
        print(f"The total cost for this configuration is ${best_config['cost']}.")
        print(f"The required delta-v for the mission (Earth orbit to Helioxis surface) is {total_required_dv} m/s.")
        print(f"This rocket's maximum delta-v is {best_config['calculated_dv']:.1f} m/s, which meets the requirement.")
        
        print("\n--- Delta-V Calculation Breakdown ---")
        # Per instructions, outputting each number in the final equation: delta_v = ve * ln(m0 / mf)
        print(f"Final Delta-V ({best_config['calculated_dv']:.1f} m/s) = Exhaust Velocity ({best_config['ve']} m/s) * ln(Total Wet Mass ({best_config['m0']:.2f} mt) / Total Dry Mass ({best_config['mf']:.2f} mt))")

        final_answer_str = f"({motor_id}) {tank_ids_str}"
        print(f"\nFinal Answer Format: {final_answer_str}")
        print(f"<<<{final_answer_str}>>>")
    else:
        print("No possible configuration meets the delta-v requirement.")
        print("<<<No solution>>>")

# Execute the function
find_cheapest_rocket()