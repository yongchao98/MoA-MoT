import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a trip from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """
    # Step 1: Determine the Required Delta-v
    # Optimal transfer delta-v (Earth to Helioxis) = 271 m/s
    # Delta-v (Helioxis surface to low orbit) = 278 m/s
    required_delta_v = 271 + 278

    # Component Data from tables
    motors = {
        1: {'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},  # ve converted to m/s
        2: {'price': 16000, 'mass': 1.54, 've': 2.3 * 1000}
    }

    tanks = {
        1: {'price': 6000, 'wet': 5.2, 'dry': 3.9},
        2: {'price': 9000, 'wet': 7.8, 'dry': 5.1},
        3: {'price': 14000, 'wet': 11.1, 'dry': 6.0},
        4: {'price': 12000, 'wet': 10.1, 'dry': 7.5}
    }

    # Step 3: Iterate Through All Possible Configurations
    min_price = float('inf')
    best_config = None
    
    tank_ids = list(tanks.keys())
    combos_3_tanks = list(itertools.combinations(tank_ids, 3))
    combos_4_tanks = list(itertools.combinations(tank_ids, 4))
    all_tank_combos = combos_3_tanks + combos_4_tanks

    # Step 4: Calculate and Compare
    for motor_id, motor_data in motors.items():
        for tank_combo in all_tank_combos:
            current_price = motor_data['price'] + sum(tanks[tid]['price'] for tid in tank_combo)

            # Calculate initial (m0) and final (mf) mass in metric tons
            m0 = motor_data['mass'] + sum(tanks[tid]['wet'] for tid in tank_combo)
            mf = motor_data['mass'] + sum(tanks[tid]['dry'] for tid in tank_combo)
            ve = motor_data['ve']
            
            # Calculate delta-v using the Tsiolkovsky rocket equation
            calculated_delta_v = ve * math.log(m0 / mf)

            # Check if the rocket meets the requirements and is cheaper
            if calculated_delta_v >= required_delta_v:
                if current_price < min_price:
                    min_price = current_price
                    best_config = {
                        'motor': motor_id,
                        'tanks': tank_combo,
                        'price': current_price,
                        'delta_v': calculated_delta_v,
                        've': ve,
                        'm0': m0,
                        'mf': mf
                    }

    # Step 5: Identify and Print the Cheapest Solution
    if best_config:
        print(f"The required delta-v for the mission is {required_delta_v} m/s.")
        print(f"The cheapest valid rocket configuration costs ${best_config['price']}.")
        
        # Output the numbers used in the final equation for the best configuration
        print("\nThe performance calculation for the chosen rocket is:")
        m0_rounded = round(best_config['m0'], 2)
        mf_rounded = round(best_config['mf'], 2)
        dv_rounded = round(best_config['delta_v'], 2)
        print(f"delta_v = {best_config['ve']} * ln({m0_rounded} / {mf_rounded}) = {dv_rounded} m/s")

        # Format the final answer string
        tank_str = ", ".join(map(str, sorted(list(best_config['tanks']))))
        final_answer_str = f"({best_config['motor']}) {tank_str};"
        
        print(f"\nThe cheapest rocket configuration is: {final_answer_str}")
        print(f"\n<<<{final_answer_str}>>>")
    else:
        print("No valid rocket configuration was found to meet the mission requirements.")

solve_rocket_problem()