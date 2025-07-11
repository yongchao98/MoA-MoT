import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a trip from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """
    # --- Data Definition ---
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 've': 1200},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 've': 2300}
    ]
    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- Step 1: Calculate Total Delta-V Requirement ---
    # Low Earth Orbit to Escape Earth (28 m/s)
    # Earth to Helioxis transfer (271 m/s)
    # Helioxis capture (low orbit to escape) (131 m/s)
    # Helioxis low orbit to surface (278 m/s)
    required_delta_v = 28 + 271 + 131 + 278

    # --- Step 3 & 4: Iterate Configurations and Find Cheapest Capable Rocket ---
    min_cost = float('inf')
    best_config_details = None

    tank_indices = list(range(len(tanks)))
    # Combinations of 3 or 4 tanks
    for num_tanks in [3, 4]:
        for tank_combo_indices in itertools.combinations(tank_indices, num_tanks):
            for motor in motors:
                current_cost = motor['price']
                m_initial = motor['wet_mass']
                m_final = motor['dry_mass']
                
                current_tanks_config = []
                for tank_idx in tank_combo_indices:
                    tank = tanks[tank_idx]
                    current_cost += tank['price']
                    m_initial += tank['wet_mass']
                    m_final += tank['dry_mass']
                    current_tanks_config.append(tank)

                # Calculate achievable delta-v (Tsiolkovsky rocket equation)
                # Ensure m_initial > m_final to avoid math domain error with log
                if m_initial > m_final:
                    achievable_delta_v = motor['ve'] * math.log(m_initial / m_final)
                else:
                    achievable_delta_v = 0
                
                # Check if it's a valid and cheaper solution
                if achievable_delta_v >= required_delta_v:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        best_config_details = {
                            'motor': motor,
                            'tanks': current_tanks_config,
                            'cost': current_cost,
                            'delta_v': achievable_delta_v,
                            'm_initial': m_initial,
                            'm_final': m_final
                        }

    # --- Step 5: Format and Display the Answer ---
    if best_config_details:
        motor_id = best_config_details['motor']['id']
        tank_ids = sorted([t['id'] for t in best_config_details['tanks']])
        result_str = f"({motor_id}) {', '.join(map(str, tank_ids))}"
        
        print(f"The total delta-v required for the mission is {required_delta_v} m/s.")
        print("\nThe cheapest capable rocket configuration is:")
        print(f"  Configuration String: {result_str}")
        print(f"  Cost: {best_config_details['cost']}")
        
        print("\nThis rocket is capable based on the Tsiolkovsky rocket equation:")
        print(f"  Δv = Ve * ln(M_initial / M_final)")
        # Here we output each number in the final equation
        print(f"  {best_config_details['delta_v']:.2f} m/s = {best_config_details['motor']['ve']} * ln({best_config_details['m_initial']:.2f} / {best_config_details['m_final']:.2f})")
        
        print(f"\n<<<{result_str}>>>")
    else:
        print("No capable rocket configuration was found within the given constraints.")
        print("\n<<<No capable rocket configuration found>>>")

solve_rocket_problem()