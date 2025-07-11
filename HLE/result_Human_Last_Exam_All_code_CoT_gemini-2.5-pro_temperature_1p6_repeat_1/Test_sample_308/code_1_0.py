import math
import itertools

def find_cheapest_rocket():
    """
    Calculates the cheapest rocket configuration for a mission from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """
    # Step 1: Define component data and mission requirements
    motors = [
        {"id": 1, "name": "Rocket motor one", "price": 8000, "mass": 1.3, "ve": 1.2 * 1000},
        {"id": 2, "name": "Rocket motor two", "price": 16000, "mass": 1.54, "ve": 2.3 * 1000}
    ]

    tanks = [
        {"id": 1, "name": "Fuel tank one", "price": 6000, "wet_mass": 5.2, "dry_mass": 3.9},
        {"id": 2, "name": "Fuel tank two", "price": 9000, "wet_mass": 7.8, "dry_mass": 5.1},
        {"id": 3, "name": "Fuel tank three", "price": 14000, "wet_mass": 11.1, "dry_mass": 6.0},
        {"id": 4, "name": "Fuel tank four", "price": 12000, "wet_mass": 10.1, "dry_mass": 7.5}
    ]

    # Delta-v from low Earth orbit to low Helioxis orbit
    delta_v_transfer = 271  # m/s
    # Delta-v from low Helioxis orbit to Helioxis surface
    delta_v_landing = 278  # m/s
    required_delta_v = delta_v_transfer + delta_v_landing

    cheapest_config = None
    min_cost = float('inf')

    # Step 2: Iterate through all possible configurations
    # One motor, with combinations of 3 or 4 tanks
    for motor in motors:
        for num_tanks in [3, 4]:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                # Step 3: Calculate cost and mass for the current configuration
                current_cost = motor["price"]
                total_wet_mass = motor["mass"]
                total_dry_mass = motor["mass"]
                tank_ids = []

                for tank in tank_combo:
                    current_cost += tank["price"]
                    total_wet_mass += tank["wet_mass"]
                    total_dry_mass += tank["dry_mass"]
                    tank_ids.append(tank["id"])
                
                # Calculate delta-v using the Tsiolkovsky rocket equation
                if total_dry_mass > 0:
                    mass_ratio = total_wet_mass / total_dry_mass
                    if mass_ratio > 1:
                        calculated_delta_v = motor["ve"] * math.log(mass_ratio)
                    else:
                        calculated_delta_v = 0
                else:
                    calculated_delta_v = 0

                # Step 4: Check if the configuration is valid and optimal
                if calculated_delta_v >= required_delta_v:
                    if current_cost < min_cost:
                        min_cost = current_cost
                        cheapest_config = {
                            "motor_id": motor["id"],
                            "tank_ids": sorted(tank_ids),
                            "cost": current_cost,
                            "delta_v": calculated_delta_v
                        }

    # Step 5: Print the final result
    if cheapest_config:
        motor_id_str = f"({cheapest_config['motor_id']})"
        tank_ids_str = ", ".join(map(str, cheapest_config['tank_ids']))
        result_string = f"{motor_id_str} {tank_ids_str}"
        print(result_string)
    else:
        print("No suitable rocket configuration was found.")

if __name__ == '__main__':
    find_cheapest_rocket()