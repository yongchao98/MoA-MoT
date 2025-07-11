import math
from itertools import combinations

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis.
    """

    # Data from tables
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

    # Step 1: Calculate required delta-v
    # low Earth orbit to low Helioxis orbit + low Helioxis orbit to Helioxis surface
    delta_v_transfer = 271  # m/s
    delta_v_landing = 278   # m/s
    required_delta_v = delta_v_transfer + delta_v_landing

    print(f"Mission requirement: Low Earth Orbit to Helioxis Surface")
    print(f"Required Delta-V = {delta_v_transfer} m/s (transfer) + {delta_v_landing} m/s (landing) = {required_delta_v} m/s\n")

    cheapest_config = None
    min_cost = float('inf')

    # Iterate through all possible numbers of tanks (3 or 4)
    for num_tanks in [3, 4]:
        # Iterate through all combinations of tanks for the given number
        for tank_combo in combinations(tanks, num_tanks):
            # Iterate through each motor
            for motor in motors:
                # Calculate properties for the current configuration
                current_cost = motor['price'] + sum(t['price'] for t in tank_combo)
                
                total_wet_mass = motor['mass'] + sum(t['wet_mass'] for t in tank_combo)
                total_dry_mass = motor['mass'] + sum(t['dry_mass'] for t in tank_combo)

                # Calculate delta-v using the Tsiolkovsky rocket equation
                # Δv = ve * ln(m_initial / m_final)
                exhaust_velocity = motor['ve']
                if total_dry_mass > 0 and total_wet_mass > total_dry_mass:
                    rocket_delta_v = exhaust_velocity * math.log(total_wet_mass / total_dry_mass)
                else:
                    rocket_delta_v = 0

                # Check if the rocket is capable and cheaper than the current best
                if rocket_delta_v >= required_delta_v:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        cheapest_config = {
                            'motor': motor,
                            'tanks': list(tank_combo),
                            'cost': current_cost,
                            'delta_v': rocket_delta_v,
                            'total_wet_mass': total_wet_mass,
                            'total_dry_mass': total_dry_mass,
                        }

    # Print the results
    if cheapest_config:
        motor = cheapest_config['motor']
        tanks_list = cheapest_config['tanks']
        tank_ids = sorted([t['id'] for t in tanks_list])

        print(f"Found cheapest capable rocket:")
        print(f"  - Configuration: Motor {motor['id']} with Tanks {', '.join(map(str, tank_ids))}")
        print(f"  - Cost: ${cheapest_config['cost']}")
        
        # Fulfilling the requirement to "output each number in the final equation"
        tank_wet_mass_str = ' + '.join(map(str, [t['wet_mass'] for t in tanks_list]))
        tank_dry_mass_str = ' + '.join(map(str, [t['dry_mass'] for t in tanks_list]))
        
        print("\nDelta-V Calculation:")
        print(f"  Delta-V = Vₑ * ln(M_initial / M_final)")
        print(f"  Delta-V = {motor['ve']} * ln(({motor['mass']} + {tank_wet_mass_str}) / ({motor['mass']} + {tank_dry_mass_str}))")
        print(f"  Delta-V = {motor['ve']} * ln({cheapest_config['total_wet_mass']:.2f} / {cheapest_config['total_dry_mass']:.2f})")
        print(f"  Delta-V = {cheapest_config['delta_v']:.2f} m/s (which is >= {required_delta_v} m/s)")
        
        final_answer_str = f"({motor['id']}) {', '.join(map(str, tank_ids))}"
        print(f"\nFinal answer format: {final_answer_str}")
        print(f"<<<{final_answer_str}>>>")
    else:
        print("No capable rocket configuration was found under the given constraints.")

solve_rocket_problem()