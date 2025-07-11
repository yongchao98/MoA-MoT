import math
from itertools import combinations

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a given space mission
    based on provided component data.
    """
    # 1. Define component data and mission requirements
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'mass': 1.54, 've': 2.3 * 1000}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet': 5.2, 'dry': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet': 7.8, 'dry': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet': 11.1, 'dry': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet': 10.1, 'dry': 7.5}
    ]

    # Required delta-v: LEO to Earth escape + Earth to Helioxis transfer + Helioxis orbit to surface
    dv_earth_escape = 28  # m/s
    dv_transfer = 271      # m/s
    dv_helioxis_landing = 278 # m/s
    dv_required = dv_earth_escape + dv_transfer + dv_helioxis_landing

    # 2. Initialize variables to track the best configuration
    min_cost = float('inf')
    best_config = None

    # 3. Iterate through all possible valid configurations
    for motor in motors:
        for num_tanks in [3, 4]:
            for tank_combo in combinations(tanks, num_tanks):
                
                # 4. Calculate total cost, wet mass, and dry mass for the current configuration
                current_cost = motor['price']
                total_wet_mass = motor['mass']
                total_dry_mass = motor['mass']

                for tank in tank_combo:
                    current_cost += tank['price']
                    total_wet_mass += tank['wet']
                    total_dry_mass += tank['dry']
                
                # 5. Calculate the rocket's delta-v using the Tsiolkovsky rocket equation
                # Ensure mass ratio is valid for the log function
                if total_dry_mass <= 0 or total_wet_mass <= total_dry_mass:
                    calculated_dv = 0
                else:
                    mass_ratio = total_wet_mass / total_dry_mass
                    calculated_dv = motor['ve'] * math.log(mass_ratio)

                # 6. Check if the rocket meets requirements and is cheaper than the current best
                if calculated_dv >= dv_required and current_cost < min_cost:
                    min_cost = current_cost
                    best_config = {
                        'motor': motor,
                        'tanks': tank_combo,
                        'cost': current_cost,
                        'dv': calculated_dv,
                        'wet_mass': total_wet_mass,
                        'dry_mass': total_dry_mass
                    }

    # 7. Print the final results
    if best_config:
        tank_ids = sorted([tank['id'] for tank in best_config['tanks']])
        tank_str = ", ".join(map(str, tank_ids))
        config_str = f"({best_config['motor']['id']}) {tank_str}"
        
        print("The cheapest rocket capable of the mission is:")
        print(f"Configuration: {config_str}")
        print(f"Total Cost: ${best_config['cost']}")
        
        print("\nThis rocket's delta-v is calculated as follows:")
        print(f"Delta-v = ve * ln(m_wet / m_dry)")
        print(f"Delta-v = {best_config['motor']['ve']} * ln({best_config['wet_mass']:.2f} / {best_config['dry_mass']:.2f})")
        print(f"Result: {best_config['dv']:.2f} m/s, which is sufficient for the required {dv_required} m/s.")
        
        print(f"<<<{config_str}>>>")
    else:
        print("No suitable rocket configuration was found to meet the mission requirements.")

# Run the function
find_cheapest_rocket()