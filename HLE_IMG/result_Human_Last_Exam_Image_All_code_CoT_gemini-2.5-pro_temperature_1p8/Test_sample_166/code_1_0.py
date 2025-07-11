import math

def calculate_reactive_power_compensation():
    """
    Calculates the required reactive power compensation from the ESS.
    """
    # Step 1: Define the loads and power factor.
    load_l1 = 1.85  # MW
    load_l2 = 1.70  # MW
    load_l3 = 1.75  # MW
    load_l4 = 1.90  # MW
    load_l5 = 2.40  # MW
    power_factor = 0.95  # lagging

    # Step 2: Calculate the total active power demand.
    total_active_power = load_l1 + load_l2 + load_l3 + load_l4 + load_l5

    # Step 3: Calculate the total reactive power demand.
    # The relationship is Q = P * tan(phi), where phi = arccos(PF).
    # This represents the total reactive power consumed by the loads.
    power_factor_angle = math.acos(power_factor)
    total_reactive_power = total_active_power * math.tan(power_factor_angle)

    # Step 4: The required compensation is equal to the total reactive power demand.
    # This is done to improve the power factor and stabilize the voltage.
    required_compensation = total_reactive_power
    
    # Print the details of the calculation
    print("--- Calculation of Required Reactive Power Compensation ---")
    print(f"Total Active Power (P_total) = {load_l1} + {load_l2} + {load_l3} + {load_l4} + {load_l5} = {total_active_power:.2f} MW")
    print(f"Load Power Factor (PF) = {power_factor}")
    print(f"Power Factor Angle (phi) = acos({power_factor}) = {math.degrees(power_factor_angle):.2f} degrees")
    print("\nTo correct the voltage drop, the ESS must compensate for the total reactive power consumed by the loads.")
    print("The required reactive power (Q_ESS) is calculated as: Q_ESS = P_total * tan(phi)")
    print(f"Q_ESS = {total_active_power:.2f} MW * tan({math.degrees(power_factor_angle):.2f}°) = {required_compensation:.2f} Mvar")
    
    # Return the final numerical answer for the platform
    # The platform expects <<<answer>>> format
    # The question does not specify the unit, but it's Mvar.
    # We will output just the number.
    print(f"\nThe required reactive power compensation from the ESS at Bus 9 is {required_compensation:.2f} Mvar.")
    

calculate_reactive_power_compensation()
print(f'<<<{3.16}>>>')
