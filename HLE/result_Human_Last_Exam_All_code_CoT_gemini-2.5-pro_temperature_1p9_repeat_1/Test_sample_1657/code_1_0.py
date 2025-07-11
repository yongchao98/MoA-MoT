from fractions import Fraction

def solve_pandora_mission():
    """
    Calculates the mission time to Pandora for both Earth and Pioneer's clocks.
    
    This function implements a physics model based on the following interpretations:
    1. Pandora's distance is derived from its recession velocity using Hubble's Law,
       assuming a Hubble constant H₀ = 70 km/s/Mpc.
    2. Pioneer's acceleration is interpreted as a percentage of Pandora's recession velocity.
    3. Motion is modeled using Newtonian kinematics, with a correction for time dilation
       to calculate ship time.
    All calculations are performed using the Fraction type for precision, simulating
    the decimal-based 'frac' type on the Wuxing computer.
    """
    
    # --- Constants ---
    # Using fractions to maintain precision as required by the Wuxing architecture simulation.
    # Speed of light in km/s
    C_LIGHT = Fraction(299792) 
    # Time conversions
    SECONDS_PER_DAY = Fraction(86400)
    SECONDS_PER_YEAR = Fraction(36525, 100) * SECONDS_PER_DAY
    # Cosmological constants for distance calculation
    HUBBLE_CONSTANT = Fraction(70)  # km/s per Mpc
    KILOMETERS_PER_MPC = Fraction(30856775814913673) # More precise value

    # --- 1. Pandora's Properties ---
    # Redshift z = (λ_obs - λ_emit) / λ_emit = (501 - 500) / 500
    redshift = Fraction(1, 500)
    # Pandora's recession velocity (v_pan = z * c)
    v_pan_kms = redshift * C_LIGHT
    # Initial distance to Pandora using Hubble's Law (D = v / H₀)
    distance_mpc = v_pan_kms / HUBBLE_CONSTANT
    d_initial_km = distance_mpc * KILOMETERS_PER_MPC
    
    # --- 2. Pioneer's Acceleration Phase (400 days) ---
    v_pioneer_kms = Fraction(40) # Initial speed
    d_accel_total_km = Fraction(0)
    t_accel_total_s = Fraction(0)
    tau_accel_total_s = Fraction(0) # Onboard proper time
    
    accel_daily_rates = [Fraction(4, 100), Fraction(2, 100), Fraction(1, 100), Fraction(5, 1000)]
    t_phase_s = Fraction(100) * SECONDS_PER_DAY

    v_start_phase = v_pioneer_kms
    for rate in accel_daily_rates:
        # Change in velocity during the 100-day phase
        # Interpretation: a = rate * v_pan [km/s/day] -> Δv = a * 100 days
        delta_v_phase = rate * v_pan_kms * 100
        v_end_phase = v_start_phase + delta_v_phase
        
        # Distance covered in phase: d = avg_v * t
        d_phase = (v_start_phase + v_end_phase) / 2 * t_phase_s
        d_accel_total_km += d_phase

        # Proper time (ship's clock) calculation for the phase
        # Δτ ≈ Δt * (1 - v_avg²/2c²)
        # More accurately: Δτ = ∫√(1-v(t)²/c²)dt ≈ ∫(1-v(t)²/2c²)dt
        # ∫v(t)²dt over [0,T] for linear v(t) is T/3 * (v_start²+v_start*v_end+v_end²)
        v_sq_integral = t_phase_s * (v_start_phase**2 + v_start_phase * v_end_phase + v_end_phase**2) / 3
        delta_tau = t_phase_s - v_sq_integral / (2 * C_LIGHT**2)
        tau_accel_total_s += delta_tau

        t_accel_total_s += t_phase_s
        v_start_phase = v_end_phase

    v_cruise_kms = v_start_phase

    # --- 3. Cruise Phase ---
    # Distance Pandora moved away during Pioneer's acceleration
    d_pan_moved_km = v_pan_kms * t_accel_total_s
    # Total distance to cover after acceleration phase
    d_rem_km = d_initial_km + d_pan_moved_km - d_accel_total_km
    # Closing speed
    v_close_kms = v_cruise_kms - v_pan_kms

    t_cruise_s = d_rem_km / v_close_kms

    # --- 4. Final Calculation ---
    # Total time as measured on Earth
    t_total_earth_s = t_accel_total_s + t_cruise_s
    total_earth_years = t_total_earth_s / SECONDS_PER_YEAR

    # Proper time (ship's clock) during cruise phase
    tau_cruise_s = t_cruise_s * (1 - (v_cruise_kms**2 / (2 * C_LIGHT**2)))
    # Total time as measured on the ship
    tau_total_ship_s = tau_accel_total_s + tau_cruise_s
    total_ship_years = tau_total_ship_s / SECONDS_PER_YEAR
    
    # --- Final Output ---
    # The C program would output integer representations of these final values.
    print(f"{int(total_earth_years)}:{int(total_ship_years)}")

solve_pandora_mission()