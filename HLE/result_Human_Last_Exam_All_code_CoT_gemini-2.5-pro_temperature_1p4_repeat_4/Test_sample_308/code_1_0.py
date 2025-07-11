import math
import itertools

def solve_rocket_problem():
    # Step 1: Define rocket component data
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'mass': 1.54, 've': 2.3 * 1000}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # Step 2: Calculate the required delta-v
    # LEO (Earth) to LEO (Helioxis) + LEO (Helioxis) to surface (Helioxis)
    dv_transfer = 271  # m/s
    dv_landing = 278   # m/s
    required_dv = dv_transfer + dv_landing
    
    # Step 3: Iterate through all valid configurations to find the cheapest one
    min_cost = float('inf')
    best_config = None

    tank_indices = list(range(len(tanks)))
    
    # Combinations of 3 or 4 tanks
    for num_tanks in [3, 4]:
        for tank_combo_indices in itertools.combinations(tank_indices, num_tanks):
            for motor in motors:
                # Calculate properties for the current configuration
                current_tanks = [tanks[i] for i in tank_combo_indices]
                
                total_price = motor['price'] + sum(t['price'] for t in current_tanks)
                
                # M_initial = motor mass + sum of wet mass of tanks
                m_initial = motor['mass'] + sum(t['wet_mass'] for t in current_tanks)
                
                # M_final = motor mass + sum of dry mass of tanks
                m_final = motor['mass'] + sum(t['dry_mass'] for t in current_tanks)

                if m_final > 0 and m_initial > m_final:
                    # Calculate achievable delta-v
                    achieved_dv = motor['ve'] * math.log(m_initial / m_final)

                    # Check if this configuration is valid and cheaper
                    if achieved_dv >= required_dv and total_price < min_cost:
                        min_cost = total_price
                        best_config = {
                            'motor': motor,
                            'tanks': current_tanks,
                            'cost': total_price,
                            'achieved_dv': achieved_dv,
                            'm_initial': m_initial,
                            'm_final': m_final,
                            'm_initial_parts': [motor['mass']] + [t['wet_mass'] for t in current_tanks],
                            'm_final_parts': [motor['mass']] + [t['dry_mass'] for t in current_tanks]
                        }

    # Step 4: Print the result
    if best_config:
        print(f"The cheapest rocket configuration is '{best_config['motor']['name']}' with tanks: {', '.join([t['name'] for t in best_config['tanks']])}.")
        print(f"Total cost: ${best_config['cost']}")
        print(f"This configuration achieves a delta-v of {best_config['achieved_dv']:.2f} m/s, which meets the requirement of {required_dv} m/s.")
        print("\nThe delta-v calculation is:")
        
        m_initial_str = ' + '.join(map(str, best_config['m_initial_parts']))
        m_final_str = ' + '.join(map(str, best_config['m_final_parts']))
        
        print(f"{best_config['achieved_dv']:.2f} = {best_config['motor']['ve']} * ln( ({m_initial_str}) / ({m_final_str}) )")
        
        tank_ids = sorted([t['id'] for t in best_config['tanks']])
        final_answer_str = f"({best_config['motor']['id']}) {', '.join(map(str, tank_ids))};"
        print(f"\n<<<({best_config['motor']['id']}) {', '.join(map(str, tank_ids))};>>>")
    else:
        print("No configuration found that can complete the mission.")

solve_rocket_problem()