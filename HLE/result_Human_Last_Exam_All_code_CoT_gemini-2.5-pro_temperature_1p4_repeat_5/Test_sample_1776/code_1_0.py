import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # Step 1: Define initial parameters
    # Half-life of Indium-111 is 2.8 days, which is 67.2 hours
    T_half_hours = 2.8 * 24

    # Desired dose activity and time of injection
    A_dose_mCi = 8.0
    t_injection = datetime(2024, 12, 23, 8, 0)  # 8:00 AM on Dec 23

    # Vial calibration details
    A_cal_mCi = 10.0
    V_cal_mL = 10.0
    t_cal = datetime(2024, 12, 26, 12, 0)  # 12:00 PM on Dec 26

    # --- Calculations ---

    # Calculate the decay constant (lambda)
    # λ = ln(2) / T_half
    decay_constant = math.log(2) / T_half_hours

    # Step 2: Calculate time difference (Δt) from injection to calibration
    time_diff = t_cal - t_injection
    time_diff_hours = time_diff.total_seconds() / 3600

    # Step 3: Calculate the vial's activity and concentration at injection time
    # Because we are calculating an activity in the past (injection) from a future
    # point (calibration), we use the inverse decay formula: A_past = A_future * e^(λ * Δt)
    activity_vial_at_injection = A_cal_mCi * math.exp(decay_constant * time_diff_hours)
    concentration_at_injection = activity_vial_at_injection / V_cal_mL

    # Step 4: Calculate the final volume to draw
    # Volume = Desired Activity / Concentration
    volume_to_draw = A_dose_mCi / concentration_at_injection

    # --- Output Results ---
    print("Step-by-step calculation of the required volume:")
    print("-" * 50)
    print(f"1. Parameters:")
    print(f"   - Indium-111 Half-life (T½): {T_half_hours:.1f} hours")
    print(f"   - Desired Dose for Patient (at 8 am): {A_dose_mCi} mCi")
    print(f"   - Vial Calibration: {A_cal_mCi} mCi in {V_cal_mL} mL at {t_cal.strftime('%I %p, %b %d')}")
    print("-" * 50)
    print(f"2. Time difference from injection (8 am, Dec 23) to calibration (12 pm, Dec 26):")
    print(f"   - Δt = {time_diff_hours:.2f} hours")
    print("-" * 50)
    print(f"3. Concentration of the vial at injection time (8 am, Dec 23):")
    print(f"   - Decay constant (λ) = ln(2) / {T_half_hours:.1f} hr = {decay_constant:.6f} hr⁻¹")
    print(f"   - Vial activity at injection = {A_cal_mCi:.1f} mCi * e^({decay_constant:.6f} * {time_diff_hours:.2f}) = {activity_vial_at_injection:.4f} mCi")
    print(f"   - Vial concentration = {activity_vial_at_injection:.4f} mCi / {V_cal_mL:.1f} mL = {concentration_at_injection:.4f} mCi/mL")
    print("-" * 50)
    print(f"4. Final Volume Calculation:")
    print(f"   The required volume is determined by the dose needed at 8 am.")
    print(f"   Final Equation: Volume = Desired Dose / Concentration at Injection")
    print(f"   Volume = {A_dose_mCi} mCi / {concentration_at_injection:.4f} mCi/mL")
    print("-" * 50)
    print(f"The volume of Indium-111 chloride that must be drawn is: {volume_to_draw:.2f} mL")


solve_radiopharmacy_calculation()
<<<3.65>>>