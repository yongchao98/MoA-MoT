import decimal

def solve_pioneer_journey():
    """
    Calculates the journey time for the Pioneer spacecraft to Pandora,
    both in Earth years and onboard (proper) years.
    """
    # Set high precision for decimal calculations to simulate the Wuxing computer's accuracy
    # and to clearly show the small difference between Earth time and onboard time.
    decimal.getcontext().prec = 50

    # --- Step 1: Define Constants and Calculate Pandora's Velocity ---
    
    # Speed of light in km/s
    c = decimal.Decimal('300000')
    
    # Calculate redshift 'z'
    lambda_obs = decimal.Decimal('501')
    lambda_rest = decimal.Decimal('500')
    z = (lambda_obs - lambda_rest) / lambda_rest
    
    # Calculate Pandora's recession velocity 'v_p' using the relativistic Doppler effect formula.
    # v/c = ((1+z)^2 - 1) / ((1+z)^2 + 1)
    one_plus_z_sq = (1 + z)**2
    v_p_over_c = (one_plus_z_sq - 1) / (one_plus_z_sq + 1)
    v_p = v_p_over_c * c

    # --- Step 2: Simulate the Journey Day by Day ---
    
    # Pioneer's initial velocity in km/s
    v_s = decimal.Decimal('40')

    # Simulation variables
    earth_days = 0
    pioneer_days = decimal.Decimal('0')

    # Pre-calculate constants for the loop
    r1 = decimal.Decimal('1.04')    # 4% acceleration
    r2 = decimal.Decimal('1.02')    # 2% acceleration
    r3 = decimal.Decimal('1.01')    # 1% acceleration
    r4 = decimal.Decimal('1.005')   # 0.5% acceleration
    one = decimal.Decimal('1')
    c_sq_2 = 2 * c**2

    # Loop until Pioneer's velocity matches/exceeds Pandora's recession velocity
    while v_s <= v_p:
        earth_days += 1
        
        # Calculate onboard time for this day using the time dilation approximation:
        # T_pioneer_day = 1_day * sqrt(1 - v_s^2/c^2) ≈ 1 - (v_s^2 / (2*c^2))
        time_adjustment = (v_s**2) / c_sq_2
        pioneer_days += (one - time_adjustment)

        # Update Pioneer's velocity for the next day based on the acceleration schedule
        if earth_days < 100:
            v_s *= r1
        elif earth_days < 200:
            v_s *= r2
        elif earth_days < 300:
            v_s *= r3
        elif earth_days < 400:
            v_s *= r4
        # After 400 days, acceleration stops and velocity becomes constant.
        
    # --- Step 3: Convert results to years and print ---
    
    days_per_year = decimal.Decimal('365.25')
    
    # a. Earth time in years
    earth_years = decimal.Decimal(earth_days) / days_per_year
    
    # b. Pioneer onboard time in years
    pioneer_years = pioneer_days / days_per_year
    
    # The final equation is the ratio of the two calculated times.
    # The problem asks for the output in the format a:b
    print(f"Pandora's recession velocity (v_p): {v_p:.10f} km/s")
    print(f"Pioneer's final velocity (v_s): {v_s:.10f} km/s")
    print(f"Total journey time in Earth days: {earth_days}")
    print(f"Total journey time in Pioneer days: {pioneer_days:.10f}")
    print("\n--- Final Answer ---")
    print("a. Time to arrive in Earth years")
    print("b. Time to arrive in onboard years")
    print(f"{earth_years:.10f}:{pioneer_years:.10f}")

solve_pioneer_journey()