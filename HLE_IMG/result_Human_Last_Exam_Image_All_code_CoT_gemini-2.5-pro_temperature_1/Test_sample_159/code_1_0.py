import math

def calculate_system_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- 1. Define Constants and Assumptions ---
    V_rf_fundamental = 1.0  # V (peak)
    f0 = 915e6  # Hz (fundamental frequency)
    R0_parasitic = 50.0  # Ohms (parasitic resistance at f0)
    R_in_rectifier = 50.0  # Ohms (assumed input impedance of the rectifier)
    voltage_drop_factor = 0.9 # Voltage drops by 10% for each higher harmonic

    harmonics = [1, 3, 5, 7]
    voltages = []
    current_voltage = V_rf_fundamental
    for _ in harmonics:
        voltages.append(current_voltage)
        current_voltage *= voltage_drop_factor

    print("--- System Parameters ---")
    print(f"Fundamental Frequency (f0): {f0/1e6} MHz")
    print(f"Rectifier Input Resistance (R_in): {R_in_rectifier} Ohms")
    print(f"Parasitic Resistance at f0 (R0): {R0_parasitic} Ohms")
    print("-" * 25)
    
    # --- 2. Calculate Power for Each Harmonic ---
    print("--- Power Calculation per Harmonic ---")
    
    p_sourced_list = []
    p_delivered_list = []

    for i, k in enumerate(harmonics):
        V_k = voltages[i]
        
        # Parasitic resistance for the k-th harmonic
        R_p_k = R0_parasitic * (k**2)
        
        # Total resistance seen by the source for this harmonic
        R_total_k = R_p_k + R_in_rectifier
        
        # Power supplied by the source for this harmonic (P = V_peak^2 / (2*R))
        p_sourced_k = (V_k**2) / (2 * R_total_k)
        p_sourced_list.append(p_sourced_k)
        
        # Power delivered to the rectifier (P = I_peak^2 * R / 2 = (V_peak/R_total)^2 * R / 2)
        p_delivered_k = (V_k**2 * R_in_rectifier) / (2 * R_total_k**2)
        p_delivered_list.append(p_delivered_k)
        
        print(f"Harmonic k={k}:")
        print(f"  Voltage (V_peak) = {V_k:.3f} V")
        print(f"  Parasitic R = {R_p_k:.1f} Ohms")
        print(f"  Total R = {R_total_k:.1f} Ohms")
        print(f"  Power Sourced = {p_sourced_k*1000:.4f} mW")
        print(f"  Power Delivered to Rectifier = {p_delivered_k*1000:.4f} mW")

    # --- 3. Calculate Total Power and Efficiency ---
    p_sourced_total = sum(p_sourced_list)
    p_delivered_total = sum(p_delivered_list)
    
    efficiency = p_delivered_total / p_sourced_total if p_sourced_total > 0 else 0

    print("-" * 25)
    print("\n--- Overall Efficiency Calculation ---")
    print("Efficiency (η) = Total Power Delivered / Total Power Sourced")
    
    # Format the numbers for the equation string
    p_del_str = " + ".join([f"{p:.4e}" for p in p_delivered_list])
    p_src_str = " + ".join([f"{p:.4e}" for p in p_sourced_list])
    
    print("\nFinal Equation:")
    print(f"η = ({p_del_str}) / ({p_src_str})")
    
    print(f"η = {p_delivered_total:.6f} W / {p_sourced_total:.6f} W")
    
    final_efficiency_percent = efficiency * 100
    print(f"η = {final_efficiency_percent:.2f}%")
    
    return final_efficiency_percent

# Execute the calculation and print the final answer in the required format
final_efficiency = calculate_system_efficiency()
print(f"\n<<<η = {final_efficiency:.2f}%>>>")