import math

def calculate_reactive_power_compensation():
    """
    Calculates the required reactive power compensation from the ESS.
    """
    # Step 1: Define the real power of each load in MW from the diagram.
    p_l1 = 1.85
    p_l2 = 1.7
    p_l3 = 1.75
    p_l4 = 1.9
    p_l5 = 2.4
    
    # List of loads for clear printing
    loads = [p_l1, p_l2, p_l3, p_l4, p_l5]
    
    # Step 2: Calculate the total real power load (P_total) by summing all loads.
    p_total = sum(loads)
    
    # Step 3: Define the system's power factor.
    power_factor = 0.95
    
    # Step 4: Calculate the total reactive power demand (Q_total).
    # The power factor angle (phi) is the arccosine of the power factor.
    # The total reactive power is P_total * tan(phi).
    phi = math.acos(power_factor)
    q_total = p_total * math.tan(phi)
    
    # The ESS must supply this amount of reactive power for compensation.
    q_ess = q_total

    print("Step 1: Calculate the total real power demand (P_total).")
    load_sum_str = " + ".join(map(str, loads))
    print(f"P_total = {load_sum_str} MW")
    print(f"P_total = {p_total:.2f} MW\n")
    
    print("Step 2: Calculate the total reactive power demand (Q_total) using the power factor (PF).")
    print(f"The formula is: Q_total = P_total * tan(acos(PF))")
    print(f"Substituting the values:")
    print(f"Q_total = {p_total:.2f} * tan(acos({power_factor}))")
    print(f"Q_total = {p_total:.2f} * {math.tan(phi):.6f}")
    print(f"Q_total = {q_total:.3f} Mvar\n")
    
    print("Step 3: The required reactive power compensation from the ESS is equal to the total reactive power demand.")
    print(f"Required Reactive Power from ESS = {q_ess:.3f} Mvar")

if __name__ == '__main__':
    calculate_reactive_power_compensation()
    
    # Calculation for the final answer
    p_total = 1.85 + 1.7 + 1.75 + 1.9 + 2.4
    power_factor = 0.95
    q_total = p_total * math.tan(math.acos(power_factor))
    print(f"\n<<<>>>") # Placeholder for the final answer wrapper
    print(f"<<<{q_total:.3f}>>>")
