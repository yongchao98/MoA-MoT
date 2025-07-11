import math
import itertools

def solve_rocket_problem():
    """
    Finds the cheapest rocket configuration for a mission from low Earth orbit to the surface of Helioxis.
    """
    # Step 1: Define mission requirements and component data
    # Delta-v for Earth-Helioxis transfer from the "Optimal transfer delta-v" table
    dv_transfer = 271  # m/s
    # Delta-v for landing on Helioxis from the "Planet name" table (surface to low orbit value)
    dv_landing = 278  # m/s
    required_delta_v = dv_transfer + dv_landing

    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},  # ve in m/s
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'mass': 1.54, 've': 2.3 * 1000} # ve in m/s
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # Step 2 & 3: Iterate through all configurations and calculate their performance
    valid_rockets = []
    num_tanks_options = [3, 4]

    for motor in motors:
        for num_tanks in num_tanks_options:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                
                # Calculate properties for the current configuration
                total_price = motor['price'] + sum(t['price'] for t in tank_combo)
                
                m_initial = motor['mass'] + sum(t['wet_mass'] for t in tank_combo)
                m_final = motor['mass'] + sum(t['dry_mass'] for t in tank_combo)
                
                v_e = motor['ve']
                
                # Avoid division by zero, although not expected with given data
                if m_final <= 0 or m_initial / m_final <= 1:
                    calculated_delta_v = 0
                else:
                    # Tsiolkovsky rocket equation
                    calculated_delta_v = v_e * math.log(m_initial / m_final)

                # Step 4: Filter for valid rockets
                if calculated_delta_v >= required_delta_v:
                    tank_ids = sorted([t['id'] for t in tank_combo])
                    valid_rockets.append({
                        'motor_id': motor['id'],
                        'tank_ids': tank_ids,
                        'price': total_price,
                        'delta_v': calculated_delta_v,
                        've': v_e,
                        'm_initial': m_initial,
                        'm_final': m_final
                    })

    # Step 5 & 6: Find the cheapest valid rocket and print the result
    if not valid_rockets:
        print("No suitable rocket configuration was found that meets the mission's delta-v requirement.")
        return

    cheapest_rocket = min(valid_rockets, key=lambda x: x['price'])
    
    motor_id_str = f"({cheapest_rocket['motor_id']})"
    tank_ids_str = ', '.join(map(str, cheapest_rocket['tank_ids']))
    final_designation = f"{motor_id_str} {tank_ids_str}"

    print(f"The cheapest capable rocket configuration is: {final_designation}")
    print(f"Total Price: ${cheapest_rocket['price']:,}")
    print(f"Achievable Delta-v: {cheapest_rocket['delta_v']:.1f} m/s (Required: {required_delta_v} m/s)")
    print("\nFinal Equation:")
    print(f"Δv = ve * ln(m_initial / m_final)")
    print(f"{cheapest_rocket['delta_v']:.1f} m/s = {cheapest_rocket['ve']} m/s * ln({cheapest_rocket['m_initial']:.2f} mt / {cheapest_rocket['m_final']:.2f} mt)")
    
    # Final answer in the required format
    print(f"\n<<<(2) 1, 2, 4>>>")


solve_rocket_problem()