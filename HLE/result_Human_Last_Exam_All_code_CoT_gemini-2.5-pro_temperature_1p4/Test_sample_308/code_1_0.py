import itertools
import math

def solve_rocket_problem():
    """
    Determines the cheapest rocket configuration (one motor and three or four fuel tanks)
    capable of traveling from low Earth orbit to the surface of Helioxis.
    """
    # 1. Define Rocket Components and Mission Parameters
    motors = {
        1: {'price': 8000, 'mass': 1.3, 've': 1.2 * 1000},  # Exhaust velocity in m/s
        2: {'price': 16000, 'mass': 1.54, 've': 2.3 * 1000}
    }

    tanks = {
        1: {'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        2: {'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        3: {'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        4: {'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    }

    # 2. Calculate Required Delta-v
    # Earth-to-Helioxis transfer + Helioxis orbit-to-surface landing
    delta_v_transfer = 271  # m/s
    delta_v_landing = 278   # m/s
    required_delta_v = delta_v_transfer + delta_v_landing
    
    # 3. Generate and Evaluate All Configurations
    cheapest_config = None
    min_cost = float('inf')

    motor_ids = list(motors.keys())
    tank_ids = list(tanks.keys())
    
    # Generate combinations of 3 and 4 tanks
    tank_combos = []
    tank_combos.extend(itertools.combinations(tank_ids, 3))
    tank_combos.extend(itertools.combinations(tank_ids, 4))

    for motor_id in motor_ids:
        for tank_combo in tank_combos:
            current_motor = motors[motor_id]

            # Calculate total cost, wet mass, and dry mass for the configuration
            current_cost = current_motor['price']
            total_wet_mass = current_motor['mass']
            total_dry_mass = current_motor['mass']

            for tank_id in tank_combo:
                current_tank = tanks[tank_id]
                current_cost += current_tank['price']
                total_wet_mass += current_tank['wet_mass']
                total_dry_mass += current_tank['dry_mass']
            
            # Use Tsiolkovsky rocket equation to calculate delta-v
            m0 = total_wet_mass  # Initial (wet) mass
            mf = total_dry_mass  # Final (dry) mass
            ve = current_motor['ve']

            if mf > 0 and m0 / mf > 1:
                calculated_delta_v = ve * math.log(m0 / mf)
            else:
                calculated_delta_v = 0

            # 4. Find the Optimal Solution
            if calculated_delta_v >= required_delta_v:
                if current_cost < min_cost:
                    min_cost = current_cost
                    cheapest_config = {
                        'motor_id': motor_id,
                        'tank_combo': sorted(list(tank_combo)),
                        'cost': current_cost,
                        'delta_v': calculated_delta_v,
                        've': ve,
                        'm0': m0,
                        'mf': mf
                    }

    # 5. Output the Result
    print(f"Required Delta-v for mission (Earth LEO to Helioxis Surface): {required_delta_v} m/s\n")
    if cheapest_config:
        motor_id_str = f"({cheapest_config['motor_id']})"
        tank_combo_str = ", ".join(map(str, cheapest_config['tank_combo']))
        final_answer = f"{motor_id_str} {tank_combo_str};"

        print("--- Cheapest Capable Rocket Configuration ---")
        print(f"Selected combination: {final_answer}")
        print(f"Total cost: ${cheapest_config['cost']}")
        print(f"Calculated Delta-v: {cheapest_config['delta_v']:.2f} m/s")
        print("\nCalculation based on the Tsiolkovsky Rocket Equation:")
        print("delta-v = Ve * ln(m0 / mf)")
        print(f"delta-v = {cheapest_config['ve']} * ln({cheapest_config['m0']:.2f} / {cheapest_config['mf']:.2f})")
    else:
        final_answer = "No suitable configuration found."
        print(final_answer)
    
    # Print the final answer for grading
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_rocket_problem()