import math

def solve_pioneer_journey():
    """
    Calculates the travel time for the Pioneer spacecraft to reach Pandora,
    both in Earth years and in onboard (Pioneer) years.
    """
    # Step 1: Define constants and calculate Pandora's properties
    # --- Physical and astronomical constants
    C_KMS = 299792.458  # Speed of light (km/s)
    S_PER_DAY = 86400.0  # Seconds per day
    DAYS_PER_YEAR = 365.25  # Days per year (accounting for leap years)
    H0_KMS_MPC = 70.0  # Hubble Constant (km/s/Mpc)
    MPC_IN_KM = 3.086e19  # Kilometers per Megaparsec

    # --- Pandora's properties from observation
    lambda_observed = 501.0  # Observed wavelength (nm)
    lambda_emitted = 500.0  # Emitted wavelength on Earth (nm)

    # Calculate redshift (z), recessional velocity, and initial distance
    z = (lambda_observed - lambda_emitted) / lambda_emitted
    v_recessional_kms = z * C_KMS
    d_initial_mpc = v_recessional_kms / H0_KMS_MPC
    d_initial_km = d_initial_mpc * MPC_IN_KM

    # Step 2: Simulate the 400-day acceleration phase
    v_kms = 40.0  # Pioneer's initial velocity (km/s)
    dist_covered_km_accel = 0.0
    onboard_time_days_accel = 0.0
    ACCEL_DURATION_DAYS = 400

    for day in range(1, ACCEL_DURATION_DAYS + 1):
        # Calculate time dilation (gamma factor) for the current day
        # gamma = 1 / sqrt(1 - v^2/c^2)
        beta_sq = (v_kms / C_KMS)**2
        gamma = 1.0 / math.sqrt(1.0 - beta_sq)

        # Update distance covered and onboard time
        dist_covered_km_accel += v_kms * S_PER_DAY
        onboard_time_days_accel += 1.0 / gamma

        # Determine acceleration multiplier for the next day's velocity
        if day <= 100:
            multiplier = 1.04
        elif day <= 200:
            multiplier = 1.02
        elif day <= 300:
            multiplier = 1.01
        else:  # day <= 400
            multiplier = 1.005
        
        v_kms *= multiplier

    # Store final values from the acceleration phase
    v_final_kms = v_kms

    # Step 3: Analyze the constant velocity phase
    # Calculate Pandora's position after the acceleration phase
    accel_duration_s = ACCEL_DURATION_DAYS * S_PER_DAY
    d_pandora_moved_km = v_recessional_kms * accel_duration_s
    
    # Calculate the remaining gap Pioneer must close
    # Gap = (Initial Distance + Pandora's Movement) - Distance Pioneer Covered
    remaining_gap_km = (d_initial_km + d_pandora_moved_km) - dist_covered_km_accel

    # Calculate the closing speed and time for this phase
    v_closing_kms = v_final_kms - v_recessional_kms
    time_const_vel_s = remaining_gap_km / v_closing_kms
    time_const_vel_days = time_const_vel_s / S_PER_DAY

    # Step 4: Calculate total times and print the final answer
    # --- Total Earth time
    total_earth_time_days = ACCEL_DURATION_DAYS + time_const_vel_days
    total_earth_time_years = total_earth_time_days / DAYS_PER_YEAR

    # --- Total Pioneer (onboard) time
    # Calculate time dilation for the constant velocity phase
    final_beta_sq = (v_final_kms / C_KMS)**2
    final_gamma = 1.0 / math.sqrt(1.0 - final_beta_sq)
    onboard_time_const_vel_days = time_const_vel_days / final_gamma
    
    total_pioneer_time_days = onboard_time_days_accel + onboard_time_const_vel_days
    total_pioneer_time_years = total_pioneer_time_days / DAYS_PER_YEAR

    # --- Output the final equations as requested
    print("Final Calculation for Earth Time:")
    print(f"Time (years) = (Acceleration Time + Constant Velocity Time) / Days per Year")
    print(f"Time (years) = ({ACCEL_DURATION_DAYS:.0f} days + {time_const_vel_days:.0f} days) / {DAYS_PER_YEAR}")
    print(f"Result: {total_earth_time_years:.0f} years")
    print("\nFinal Calculation for Pioneer Onboard Time:")
    print(f"Time (years) = (Onboard Accel Time + Onboard Const-Vel Time) / Days per Year")
    print(f"Time (years) = ({onboard_time_days_accel:.0f} days + {onboard_time_const_vel_days:.0f} days) / {DAYS_PER_YEAR}")
    print(f"Result: {total_pioneer_time_years:.0f} years")
    
    # Final answer in the required a:b format
    a = round(total_earth_time_years)
    b = round(total_pioneer_time_years)
    
    print(f"\n<<<Answer>>>\n{a}:{b}")


solve_pioneer_journey()