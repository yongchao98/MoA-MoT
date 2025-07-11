import math
import itertools

def solve_rocket_problem():
    """
    Calculates the cheapest rocket configuration for a trip from low Earth orbit
    to the surface of Helioxis, given a set of components and constraints.
    """
    # --- Data Definition ---
    motors = [
        {'id': 1, 'name': 'Rocket motor one', 'price': 8000, 'wet_mass': 1.3, 'dry_mass': 1.3, 've': 1.2},
        {'id': 2, 'name': 'Rocket motor two', 'price': 16000, 'wet_mass': 1.54, 'dry_mass': 1.54, 've': 2.3}
    ]

    tanks = [
        {'id': 1, 'name': 'Fuel tank one', 'price': 6000, 'wet_mass': 5.2, 'dry_mass': 3.9},
        {'id': 2, 'name': 'Fuel tank two', 'price': 9000, 'wet_mass': 7.8, 'dry_mass': 5.1},
        {'id': 3, 'name': 'Fuel tank three', 'price': 14000, 'wet_mass': 11.1, 'dry_mass': 6.0},
        {'id': 4, 'name': 'Fuel tank four', 'price': 12000, 'wet_mass': 10.1, 'dry_mass': 7.5}
    ]

    # --- Step 1: Calculate Required Delta-v ---
    # Optimal transfer delta-v (Earth to Helioxis) + Delta-v (Helioxis orbit to surface)
    dv_transfer = 271  # m/s
    dv_landing = 278   # m/s
    required_dv = dv_transfer + dv_landing

    # --- Step 2-4: Iterate through configurations and find the cheapest viable option ---
    min_cost = float('inf')
    best_config_details = None

    # Constraint: 1 motor, and 3 or 4 fuel tanks
    for motor in motors:
        for num_tanks in [3, 4]:
            for tank_combo in itertools.combinations(tanks, num_tanks):
                
                # a. Calculate total cost
                current_cost = motor['price'] + sum(t['price'] for t in tank_combo)

                # Only proceed if it's potentially cheaper than the best found so far
                if current_cost < min_cost:
                    # b. Calculate total mass
                    total_wet_mass = motor['wet_mass'] + sum(t['wet_mass'] for t in tank_combo)
                    total_dry_mass = motor['dry_mass'] + sum(t['dry_mass'] for t in tank_combo)
                    
                    # c. Calculate achievable delta-v
                    exhaust_velocity_ms = motor['ve'] * 1000
                    mass_ratio = total_wet_mass / total_dry_mass
                    calculated_dv = exhaust_velocity_ms * math.log(mass_ratio)
                    
                    # d. Check if the rocket is viable
                    if calculated_dv >= required_dv:
                        min_cost = current_cost
                        best_config_details = {
                            "motor": motor,
                            "tanks": tank_combo,
                            "cost": current_cost,
                            "wet_mass": total_wet_mass,
                            "dry_mass": total_dry_mass,
                            "dv": calculated_dv
                        }

    # --- Step 5: Format and Print the Final Answer ---
    if best_config_details:
        b_motor = best_config_details["motor"]
        b_tanks = best_config_details["tanks"]
        
        # Required Delta-v breakdown
        print(f"Journey: Low Earth Orbit to Helioxis Surface")
        print(f"Required Delta-v = (Earth-Helioxis Transfer) + (Helioxis Landing)")
        print(f"Required Delta-v = {dv_transfer} m/s + {dv_landing} m/s = {required_dv} m/s\n")
        
        # Best configuration found
        tank_ids = sorted([str(t['id']) for t in b_tanks])
        print(f"Best Configuration Found: Motor {b_motor['id']} with Tanks {', '.join(tank_ids)}")

        # Cost Equation
        tank_prices = [t['price'] for t in b_tanks]
        tank_price_str = " + ".join(map(str, tank_prices))
        print(f"Cost = {b_motor['price']} (Motor) + {tank_price_str} (Tanks) = {best_config_details['cost']}")

        # Delta-v Equation
        print(f"Achievable Delta-v = Ve * ln(M0 / Mf)")
        print(f"Achievable Delta-v = {b_motor['ve'] * 1000} m/s * ln({best_config_details['wet_mass']:.2f} / {best_config_details['dry_mass']:.2f}) = {best_config_details['dv']:.2f} m/s")
        
        # Final Answer
        final_answer_str = f"({b_motor['id']}) {', '.join(tank_ids)}"
        print(f"\nThe cheapest rocket capable of the journey is: {final_answer_str}")
        print(f"\n<<<{final_answer_str}>>>")
    else:
        print("No suitable rocket configuration could be found to meet the delta-v requirement.")

solve_rocket_problem()