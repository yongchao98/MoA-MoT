import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a given space mission
    and prints the result.
    """
    # --- Data from the tables ---
    motors = [
        {'id': 1, 'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},  # ve converted to m/s
        {'id': 2, 'price': 16000, 'mass': 1.54, 've': 2.3 * 1000}, # ve converted to m/s
    ]

    tanks = [
        {'id': 1, 'price': 6000, 'wet': 5.2, 'dry': 3.9},
        {'id': 2, 'price': 9000, 'wet': 7.8, 'dry': 5.1},
        {'id': 3, 'price': 14000, 'wet': 11.1, 'dry': 6.0},
        {'id': 4, 'price': 12000, 'wet': 10.1, 'dry': 7.5},
    ]

    # --- Step 1: Calculate total required Delta-V ---
    dv_earth_escape = 28  # low orbit to escape for Earth
    dv_transfer_earth_helioxis = 271 # Earth to Helioxis
    dv_helioxis_landing = 278 # surface to low orbit for Helioxis
    
    required_dv = dv_earth_escape + dv_transfer_earth_helioxis + dv_helioxis_landing
    
    # --- Step 2: Iterate through all configurations ---
    min_cost = float('inf')
    best_config = None
    best_config_details = {}

    # Get combinations of 3 and 4 tanks
    tank_combinations = []
    tank_combinations.extend(itertools.combinations(tanks, 3))
    tank_combinations.extend(itertools.combinations(tanks, 4))

    for motor in motors:
        for tank_combo in tank_combinations:
            # Calculate properties of the current configuration
            current_cost = motor['price'] + sum(t['price'] for t in tank_combo)
            
            # Optimization: if already more expensive, skip
            if current_cost >= min_cost:
                continue
                
            total_tank_wet_mass = sum(t['wet'] for t in tank_combo)
            total_tank_dry_mass = sum(t['dry'] for t in tank_combo)

            # m_0 = initial mass, m_f = final mass
            m0 = motor['mass'] + total_tank_wet_mass
            mf = motor['mass'] + total_tank_dry_mass
            
            # Calculate achievable Delta-V using Tsiolkovsky rocket equation
            if mf <= 0 or m0 / mf <= 1:
                achieved_dv = 0
            else:
                achieved_dv = motor['ve'] * math.log(m0 / mf)
            
            # Check if this configuration is sufficient and cheaper
            if achieved_dv >= required_dv and current_cost < min_cost:
                min_cost = current_cost
                best_config = (motor, tank_combo)
                best_config_details = {
                    'm0': m0,
                    'mf': mf,
                    've': motor['ve'],
                    'cost': current_cost,
                    'dv': achieved_dv
                }

    # --- Step 3: Print the results ---
    if best_config:
        motor, tank_combo = best_config
        tank_ids_str = ", ".join(sorted([str(t['id']) for t in tank_combo]))
        
        print(f"Mission Requirement: Travel from low Earth orbit to Helioxis surface.")
        print(f"Required Delta-v = {dv_earth_escape} (Earth escape) + {dv_transfer_earth_helioxis} (transfer) + {dv_helioxis_landing} (Helioxis landing) = {required_dv} m/s\n")
        print(f"Cheapest sufficient configuration found:")
        print(f"  Rocket: Motor {motor['id']} with tanks {tank_ids_str}")
        print(f"  Total Cost: ${min_cost}\n")
        print(f"Verification using the Tsiolkovsky Rocket Equation:")
        print(f"  Δv = ve * ln(m0 / mf)")
        print(f"  Numbers for the final equation:")
        print(f"    Exhaust velocity (ve): {best_config_details['ve']} m/s")
        print(f"    Initial wet mass (m0): {motor['mass']} (motor) + {round(best_config_details['m0'] - motor['mass'], 2)} (tanks) = {round(best_config_details['m0'], 2)} mt")
        print(f"    Final dry mass (mf): {motor['mass']} (motor) + {round(best_config_details['mf'] - motor['mass'], 2)} (tanks) = {round(best_config_details['mf'], 2)} mt")
        print(f"  Calculation: Δv = {best_config_details['ve']} * ln({round(best_config_details['m0'], 2)} / {round(best_config_details['mf'], 2)}) = {round(best_config_details['dv'])} m/s")
        print(f"  Result: {round(best_config_details['dv'])} m/s is greater than or equal to the required {required_dv} m/s.\n")
        print("Final Answer:")
        print(f"({motor['id']}) {tank_ids_str}")

    else:
        print("No suitable rocket configuration found within the given constraints.")

# Run the analysis
find_cheapest_rocket()