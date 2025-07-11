import math

def solve_pandora_mission():
    """
    Calculates the travel time to Pandora for the Pioneer probe.
    """
    # --- Constants ---
    C_KMS = 299792.458  # Speed of light in km/s
    V0_KMS = 40.0         # Initial speed in km/s
    SECS_PER_DAY = 86400.0
    DAYS_PER_YEAR = 365.25
    ACCEL_DURATION_DAYS = 400

    # --- Physics Calculations ---
    # Pandora's recessional velocity based on redshift z = (501-500)/500 = 1/500
    v_rec_kms = C_KMS / 500.0

    # --- Part 1: Simulate the 400-day acceleration phase ---
    v_kms = V0_KMS
    dist_covered_accel_km = 0.0
    time_onboard_accel_s = 0.0

    for day in range(ACCEL_DURATION_DAYS):
        # Determine the acceleration factor for the current day
        if day < 100:
            k = 1.04
        elif day < 200:
            k = 1.02
        elif day < 300:
            k = 1.01
        else:  # day < 400
            k = 1.005

        # Calculate time dilation (gamma) for the current velocity
        # Use velocity at the start of the day as an approximation for the whole day
        beta_sq = (v_kms**2) / (C_KMS**2)
        gamma = 1.0 / math.sqrt(1.0 - beta_sq)

        # Accumulate onboard time (dilated time)
        time_onboard_accel_s += SECS_PER_DAY / gamma

        # Accumulate distance covered
        dist_covered_accel_km += v_kms * SECS_PER_DAY

        # Update velocity for the next day
        v_kms *= k

    v_final_kms = v_kms

    # --- Part 2: Calculate total Earth time ---
    # Based on the assumption that the initial distance is the distance covered during acceleration,
    # the total time T is given by the chase equation: T = (v_final * t_accel) / (v_final - v_rec)
    
    print("Calculating total Earth time using the equation: T = (v_final * t_accel) / (v_final - v_rec)")
    print(f"v_final = {v_final_kms:.2f} km/s")
    print(f"t_accel = {ACCEL_DURATION_DAYS} days")
    print(f"v_rec = {v_rec_kms:.2f} km/s")
    
    if v_final_kms <= v_rec_kms:
        print("Error: Final velocity is not sufficient to catch up to Pandora.")
        return

    total_earth_time_days = (v_final_kms * ACCEL_DURATION_DAYS) / (v_final_kms - v_rec_kms)
    total_earth_time_years = total_earth_time_days / DAYS_PER_YEAR

    # --- Part 3: Calculate total onboard time ---
    # Duration of the constant-velocity chase phase, in Earth time
    chase_phase_earth_days = total_earth_time_days - ACCEL_DURATION_DAYS
    chase_phase_earth_s = chase_phase_earth_days * SECS_PER_DAY

    # Onboard time for the chase phase (dilated)
    final_beta_sq = (v_final_kms**2) / (C_KMS**2)
    final_gamma = 1.0 / math.sqrt(1.0 - final_beta_sq)
    chase_phase_onboard_s = chase_phase_earth_s / final_gamma

    # Total onboard time is the sum of the two phases
    total_onboard_time_s = time_onboard_accel_s + chase_phase_onboard_s
    total_onboard_time_years = total_onboard_time_s / (SECS_PER_DAY * DAYS_PER_YEAR)

    # --- Final Answer ---
    # Rounding to the nearest integer as is common for such problems
    answer_a = int(round(total_earth_time_years))
    answer_b = int(round(total_onboard_time_years))

    print("\n--- Results ---")
    print(f"a. Earth time to arrival: {answer_a} years")
    print(f"b. Onboard time to arrival: {answer_b} years")


solve_pandora_mission()