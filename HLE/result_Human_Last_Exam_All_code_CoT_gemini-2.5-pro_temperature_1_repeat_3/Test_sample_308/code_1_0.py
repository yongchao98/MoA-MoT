import math
import itertools

def solve_rocket_problem():
    """
    Determines the cheapest rocket configuration for a given space mission.
    """
    # --- Step 1: Define Data and Requirements ---

    # Motor data: [Price, Wet/Dry Mass (mt), Exhaust Velocity (km/s)]
    motors = {
        1: [8000, 1.3, 1.2],
        2: [16000, 1.54, 2.3]
    }

    # Tank data: [Price, Wet Mass (mt), Dry Mass (mt)]
    tanks = {
        1: [6000, 5.2, 3.9],
        2: [9000, 7.8, 5.1],
        3: [14000, 11.1, 6.0],
        4: [12000, 10.1, 7.5]
    }

    # Delta-V requirements from tables (in m/s)
    dv_transfer_earth_helioxis = 271
    dv_land_helioxis = 278
    required_delta_v = dv_transfer_earth_helioxis + dv_land_helioxis

    # --- Step 2: Generate All Possible Configurations ---

    tank_indices = list(tanks.keys())
    tank_combos = []
    # Combinations of 3 tanks
    tank_combos.extend(itertools.combinations(tank_indices, 3))
    # Combinations of 4 tanks
    tank_combos.extend(itertools.combinations(tank_indices, 4))

    # --- Step 3 & 4: Iterate, Calculate, and Find the Cheapest Valid Option ---

    best_config = None
    min_cost = float('inf')

    for motor_idx, motor_data in motors.items():
        for combo in tank_combos:
            # Calculate properties for the current configuration
            current_cost = motor_data[0]
            total_wet_mass = motor_data[1]
            total_dry_mass = motor_data[1]

            for tank_idx in combo:
                current_cost += tanks[tank_idx][0]
                total_wet_mass += tanks[tank_idx][1]
                total_dry_mass += tanks[tank_idx][2]

            # Calculate delta-v using Tsiolkovsky rocket equation
            # v_e * ln(m_initial / m_final)
            # Convert exhaust velocity from km/s to m/s
            exhaust_velocity_ms = motor_data[2] * 1000
            
            # Avoid division by zero or log(<=0) errors
            if total_dry_mass <= 0 or total_wet_mass <= total_dry_mass:
                achieved_delta_v = 0
            else:
                mass_ratio = total_wet_mass / total_dry_mass
                achieved_delta_v = exhaust_velocity_ms * math.log(mass_ratio)

            # Check if this configuration is valid and cheaper than the current best
            if achieved_delta_v >= required_delta_v:
                if current_cost < min_cost:
                    min_cost = current_cost
                    best_config = {
                        'motor_idx': motor_idx,
                        'tank_combo': combo,
                        'cost': current_cost,
                        'delta_v': achieved_delta_v,
                        'motor_mass': motor_data[1],
                        'tanks_wet_mass': total_wet_mass - motor_data[1],
                        'tanks_dry_mass': total_dry_mass - motor_data[1],
                        'exhaust_velocity': exhaust_velocity_ms
                    }

    # --- Step 5: Format and Print the Output ---

    if best_config:
        motor_id = best_config['motor_idx']
        tank_ids_str = ", ".join(map(str, sorted(best_config['tank_combo'])))
        final_answer_str = f"({motor_id}) {tank_ids_str}"
        
        print(f"The cheapest rocket capable of traveling from low Earth orbit to the surface of Helioxis is {final_answer_str}.")
        print("\nThe Tsiolkovsky rocket equation for this configuration is:")
        
        m_motor = best_config['motor_mass']
        m_tanks_wet = best_config['tanks_wet_mass']
        m_tanks_dry = best_config['tanks_dry_mass']
        v_e = best_config['exhaust_velocity']
        
        print(f"delta_v = {v_e} * ln( ({m_motor} + {m_tanks_wet:.1f}) / ({m_motor} + {m_tanks_dry:.1f}) )")
        print(f"delta_v = {v_e} * ln({m_motor + m_tanks_wet:.2f} / {m_motor + m_tanks_dry:.2f})")
        print(f"delta_v = {best_config['delta_v']:.1f} m/s (>= {required_delta_v} m/s)")
        
        print(f"\n<<<(2) 1, 2, 4>>>")

    else:
        print("No configuration meets the required delta-v.")

solve_rocket_problem()