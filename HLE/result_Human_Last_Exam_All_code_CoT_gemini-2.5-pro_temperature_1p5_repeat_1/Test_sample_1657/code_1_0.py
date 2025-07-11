import math

def calculate_travel_times():
    """
    Calculates the travel time to Pandora for Pioneer from both Earth's 
    and the onboard clock's perspectives.
    """
    # --- Constants and Initial Parameters ---
    # The speed of light in km/s
    C_KMS = 299792.458 
    # Time conversions
    SECONDS_PER_DAY = 86400.0
    DAYS_PER_YEAR = 365.25
    SECONDS_PER_YEAR = SECONDS_PER_DAY * DAYS_PER_YEAR
    # Distance conversion: km per light-year
    KM_PER_LY = C_KMS * SECONDS_PER_YEAR

    # --- Key Assumption ---
    # The distance to Pandora is assumed to be 500 light-years, based on the numbers 
    # in the problem statement (500 nm), as it's otherwise not specified.
    DISTANCE_LY = 500.0
    DISTANCE_KM = DISTANCE_LY * KM_PER_LY

    # Pioneer's initial parameters
    INITIAL_V_KMS = 40.0
    ACCEL_DURATION_DAYS = 400

    # --- Phase 1: Acceleration Simulation (First 400 days) ---
    v_kms = INITIAL_V_KMS
    d_accel_km = 0.0
    t_pioneer_accel_s = 0.0

    for day in range(1, ACCEL_DURATION_DAYS + 1):
        v_start_day = v_kms
        
        # Determine the velocity multiplier based on the day
        if 1 <= day <= 100:
            multiplier = 1.04
        elif 101 <= day <= 200:
            multiplier = 1.02
        elif 201 <= day <= 300:
            multiplier = 1.01
        else:  # 301 to 400
            multiplier = 1.005
            
        v_end_day = v_start_day * multiplier
        v_avg_day = (v_start_day + v_end_day) / 2.0
        
        # Add distance covered this day
        d_accel_km += v_avg_day * SECONDS_PER_DAY
        
        # Calculate time dilation (Lorentz factor gamma) for the day
        # and add to Pioneer's elapsed time (proper time).
        beta_sq_avg = (v_avg_day / C_KMS)**2
        gamma_avg = 1.0 / math.sqrt(1.0 - beta_sq_avg)
        t_pioneer_accel_s += SECONDS_PER_DAY / gamma_avg
        
        v_kms = v_end_day

    v_final_kms = v_kms
    t_earth_accel_s = float(ACCEL_DURATION_DAYS) * SECONDS_PER_DAY

    # --- Phase 2: Coasting at Constant Velocity ---
    d_rem_km = DISTANCE_KM - d_accel_km
    t_coast_earth_s = d_rem_km / v_final_kms
    
    # Calculate time dilation for the coasting phase
    beta_final_sq = (v_final_kms / C_KMS)**2
    gamma_final = 1.0 / math.sqrt(1.0 - beta_final_sq)
    t_coast_pioneer_s = t_coast_earth_s / gamma_final

    # --- Phase 3: Total Time Calculation and Final Output ---
    total_earth_time_s = t_earth_accel_s + t_coast_earth_s
    total_pioneer_time_s = t_pioneer_accel_s + t_coast_pioneer_s
    
    total_earth_years = total_earth_time_s / SECONDS_PER_YEAR
    total_pioneer_years = total_pioneer_time_s / SECONDS_PER_YEAR

    # Outputting the numbers in the final equations as requested
    print("--- Final Equation Components ---")
    print(f"Earth Time (years) = (Acceleration Time + Coasting Time) / Seconds per Year")
    print(f"Earth Time (years) = ({t_earth_accel_s:.2f} s + {t_coast_earth_s:.2f} s) / {SECONDS_PER_YEAR:.2f} s/yr\n")
    
    print(f"Pioneer Time (years) = (Acceleration Time + Coasting Time) / Seconds per Year")
    print(f"Pioneer Time (years) = ({t_pioneer_accel_s:.2f} s + {t_coast_pioneer_s:.2f} s) / {SECONDS_PER_YEAR:.2f} s/yr\n")
    
    # Final calculated answer
    answer_a = round(total_earth_years)
    answer_b = round(total_pioneer_years)
    
    print("--- Final Answers ---")
    print(f"a. Years to arrive (Earth time): {answer_a}")
    print(f"b. Years to arrive (Pioneer time): {answer_b}")

if __name__ == '__main__':
    calculate_travel_times()