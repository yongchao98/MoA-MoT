import math
import itertools

def solve_rocket_problem():
    """
    Finds the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis.
    """
    # 1. Define mission requirements and component data
    
    # Delta-v from low Earth orbit to low Helioxis orbit
    dv_transfer = 271  # m/s
    # Delta-v from low Helioxis orbit to its surface
    dv_landing = 278  # m/s
    REQUIRED_DV = dv_transfer + dv_landing

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

    # 2. Initialize variables to store the best configuration found
    cheapest_config = {
        'price': float('inf'),
        'motor': None,
        'tanks': None,
        'delta_v': 0,
        'm_initial': 0,
        'm_final': 0
    }

    # 3. Iterate through all possible valid configurations
    tank_indices = list(range(len(tanks)))
    
    for motor in motors:
        # Combinations of 3 or 4 tanks
        for num_tanks in [3, 4]:
            for tank_combo_indices in itertools.combinations(tank_indices, num_tanks):
                
                # a. Calculate total price and mass for the current configuration
                current_tanks = [tanks[i] for i in tank_combo_indices]
                
                total_price = motor['price'] + sum(t['price'] for t in current_tanks)

                # Optimization: if it's already more expensive, skip calculation
                if total_price >= cheapest_config['price']:
                    continue
                    
                m_initial = motor['mass'] + sum(t['wet_mass'] for t in current_tanks)
                m_final = motor['mass'] + sum(t['dry_mass'] for t in current_tanks)

                # b. Calculate delta-v using the Tsiolkovsky rocket equation
                if m_final > 0 and m_initial > m_final:
                    # Mass ratio must be > 1
                    mass_ratio = m_initial / m_final
                    calculated_dv = motor['ve'] * math.log(mass_ratio)
                else:
                    calculated_dv = 0

                # c. Check if the rocket meets requirements and is the cheapest so far
                if calculated_dv >= REQUIRED_DV:
                    cheapest_config = {
                        'price': total_price,
                        'motor': motor,
                        'tanks': current_tanks,
                        'delta_v': calculated_dv,
                        'm_initial': m_initial,
                        'm_final': m_final
                    }

    # 4. Print the final result
    if cheapest_config['motor']:
        best_motor = cheapest_config['motor']
        best_tanks = cheapest_config['tanks']
        tank_ids_str = ", ".join(str(t['id']) for t in sorted(best_tanks, key=lambda x: x['id']))
        config_str = f"({best_motor['id']}) {tank_ids_str}"
        
        print("Cheapest capable rocket configuration found:")
        print(f"Configuration: {config_str}")
        print(f"Total Price: ${cheapest_config['price']}")
        print(f"Required Delta-v: {REQUIRED_DV} m/s")
        print(f"Achieved Delta-v: {cheapest_config['delta_v']:.2f} m/s")
        print("\nDelta-v Calculation:")
        print(f"Δv = Ve * ln(M_initial / M_final)")
        print(f"Δv = {best_motor['ve']} * ln({cheapest_config['m_initial']:.2f} / {cheapest_config['m_final']:.2f})")

        # The final answer in the required format
        print(f"\n<<<{config_str}>>>")
    else:
        print("No configuration found that meets the mission requirements.")

solve_rocket_problem()