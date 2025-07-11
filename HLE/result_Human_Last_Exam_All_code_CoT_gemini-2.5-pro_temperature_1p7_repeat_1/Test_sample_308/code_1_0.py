import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis, based on the provided component and travel data.
    """
    # --- Step 1: Define Data ---

    # Motor data: id, price, wet mass (mt), dry mass (mt), exhaust velocity (m/s)
    motors = [
        {'id': 1, 'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 'exhaust_velocity': 1200},
        {'id': 2, 'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 'exhaust_velocity': 2300}
    ]

    # Tank data: id, price, wet mass (mt), dry mass (mt)
    tanks = [
        {'id': 1, 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # Mission delta-v requirements (m/s)
    delta_v_transfer_earth_to_helioxis = 271
    delta_v_land_on_helioxis = 278

    # --- Step 2: Calculate Total Required Delta-v ---
    
    required_delta_v = delta_v_transfer_earth_to_helioxis + delta_v_land_on_helioxis
    print(f"Calculating required delta-v for the mission:")
    print(f"Delta-v (Earth to Helioxis transfer) + Delta-v (Helioxis landing) = Required Delta-v")
    print(f"{delta_v_transfer_earth_to_helioxis} m/s + {delta_v_land_on_helioxis} m/s = {required_delta_v} m/s\n")

    # --- Step 3: Iterate Through All Configurations to Find the Cheapest Valid One ---

    min_cost = float('inf')
    best_config = None
    best_config_details = {}

    # Generate combinations of 3 and 4 tanks
    tank_combos_3 = list(itertools.combinations(tanks, 3))
    tank_combos_4 = list(itertools.combinations(tanks, 4))
    all_tank_combos = tank_combos_3 + tank_combos_4

    for motor in motors:
        for tank_combo in all_tank_combos:
            # Calculate total cost
            tank_cost = sum(tank['price'] for tank in tank_combo)
            total_cost = motor['price'] + tank_cost

            # Optimization: if it's already more expensive than our best find, skip it
            if total_cost >= min_cost:
                continue

            # Calculate total mass
            tank_wet_mass = sum(tank['wet_mass'] for tank in tank_combo)
            tank_dry_mass = sum(tank['dry_mass'] for tank in tank_combo)

            m_initial = motor['wet_mass'] + tank_wet_mass
            m_final = motor['dry_mass'] + tank_dry_mass

            # Calculate achievable delta-v using the Tsiolkovsky rocket equation
            # Δv = v_e * ln(m_initial / m_final)
            if m_initial > m_final:
                 achieved_delta_v = motor['exhaust_velocity'] * math.log(m_initial / m_final)
            else:
                 achieved_delta_v = 0
            
            # Check if this configuration is capable and cheaper
            if achieved_delta_v >= required_delta_v:
                if total_cost < min_cost:
                    min_cost = total_cost
                    best_config_details = {
                        'motor_id': motor['id'],
                        'tank_ids': sorted([tank['id'] for tank in tank_combo]),
                        'cost': total_cost,
                        'delta_v': achieved_delta_v,
                        'v_e': motor['exhaust_velocity'],
                        'm_initial': m_initial,
                        'm_final': m_final
                    }

    # --- Step 4: Display the Result ---

    if best_config_details:
        print("Found the cheapest capable rocket configuration:")
        
        # Display the equation for the winning configuration
        b = best_config_details
        print(f"Tsiolkovsky equation for the best configuration:")
        print(f"Δv = exhaust_velocity * ln(initial_mass / final_mass)")
        print(f"Achieved Δv = {b['v_e']} * ln({b['m_initial']:.2f} / {b['m_final']:.2f}) = {b['delta_v']:.2f} m/s\n")
        
        # Format the final answer as requested
        motor_id_str = f"({b['motor_id']})"
        tank_ids_str = ', '.join(map(str, b['tank_ids']))
        final_answer = f"{motor_id_str} {tank_ids_str}"
        
        print(f"The cheapest capable configuration is: {final_answer}")
        print(f"Total cost: ${b['cost']}")
        
        print(f"\n<<<{final_answer}>>>")
    else:
        print("No capable rocket configuration was found within the given constraints.")

# Run the analysis
find_cheapest_rocket()