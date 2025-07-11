import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride to be drawn for a patient dose.
    """
    # --- 1. Define constants and initial values ---
    A_req_inj = 8.0     # Required activity at injection time (mCi)
    A_stock_cal = 10.0    # Stock activity at calibration time (mCi)
    V_stock = 10.0      # Stock volume (mL)
    T_half_days = 2.80  # Half-life of Indium-111 in days

    # --- 2. Define time points ---
    # The year 2024 is specified in the problem.
    t_inj = datetime(2024, 12, 23, 8, 0)   # Injection time: Dec 23, 8:00 am
    t_cal = datetime(2024, 12, 26, 12, 0)  # Calibration time: Dec 26, 12:00 pm

    # --- 3. Calculate time difference in hours ---
    time_difference = t_cal - t_inj
    delta_t_hours = time_difference.total_seconds() / 3600.0

    # --- 4. Calculate decay constant (lambda) in hours^-1 ---
    T_half_hours = T_half_days * 24.0
    lambda_val = math.log(2) / T_half_hours

    # --- 5. Calculate stock concentration at calibration time ---
    C_stock_cal = A_stock_cal / V_stock

    # --- 6. Calculate the required activity at the calibration time ---
    # A_req_cal = A_req_inj * e^(-lambda * t)
    A_req_cal = A_req_inj * math.exp(-lambda_val * delta_t_hours)

    # --- 7. Calculate the volume to draw ---
    V_draw = A_req_cal / C_stock_cal

    # --- 8. Print the calculation steps and final answer ---
    print("--- Problem Setup ---")
    print(f"Required Dose Activity at Injection: {A_req_inj:.1f} mCi")
    print(f"Stock Vial Calibration: {A_stock_cal:.1f} mCi in {V_stock:.1f} mL")
    print(f"Indium-111 Half-Life: {T_half_days} days ({T_half_hours:.1f} hours)")
    
    print("\n--- Calculation Steps ---")
    print(f"1. Time from Injection to Calibration (t):")
    print(f"   t = (Time at Calibration) - (Time at Injection)")
    print(f"   t = ({t_cal.strftime('%Y-%m-%d %H:%M')}) - ({t_inj.strftime('%Y-%m-%d %H:%M')})")
    print(f"   t = {delta_t_hours:.2f} hours")
    
    print(f"\n2. Decay Constant (λ):")
    print(f"   λ = ln(2) / Half-Life")
    print(f"   λ = {math.log(2):.6f} / {T_half_hours:.1f} hours = {lambda_val:.6f} hr⁻¹")

    print(f"\n3. Dose Activity Needed at Calibration Time (A_needed_cal):")
    print(f"   A_needed_cal = (Dose at Injection) * e^(-λ*t)")
    print(f"   A_needed_cal = {A_req_inj:.1f} mCi * e^(-{lambda_val:.6f} * {delta_t_hours:.2f})")
    print(f"   A_needed_cal = {A_req_cal:.4f} mCi")

    print(f"\n4. Stock Concentration at Calibration Time (C_stock_cal):")
    print(f"   C_stock_cal = (Stock Activity at Cal) / (Stock Volume)")
    print(f"   C_stock_cal = {A_stock_cal:.1f} mCi / {V_stock:.1f} mL = {C_stock_cal:.1f} mCi/mL")

    print("\n--- Final Volume Calculation ---")
    print("Volume to Draw = (Dose Activity Needed at Cal) / (Stock Concentration at Cal)")
    print(f"Volume to Draw = {A_req_cal:.4f} mCi / {C_stock_cal:.1f} mCi/mL")
    print(f"Result: The volume of Indium-111 chloride to be drawn is {V_draw:.2f} mL.")
    
    print("\n--- Final Equation with all numbers plugged in ---")
    # Equation: V = (A_req_inj * e^(-(ln(2)/T_half_hr) * t_hr)) / (A_stock_cal / V_stock)
    final_equation = f"Volume (mL) = ({A_req_inj:.1f} * exp(-({math.log(2):.4f} / {T_half_hours:.1f}) * {delta_t_hours:.2f})) / ({A_stock_cal:.1f} / {V_stock:.1f})"
    print(final_equation)
    return V_draw

# Run the calculation and store the result
final_volume = solve_radiopharmacy_calculation()
# The final answer tag will be appended after the code block by the system.
# The value is final_volume, which is approximately 3.65