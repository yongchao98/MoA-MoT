import math
from fractions import Fraction

def solve_pandora_mission():
    """
    This program simulates the Pioneer mission to Pandora, respecting the
    computational constraints of the Wuxing architecture by using the Fraction
    type to avoid floating-point numbers. It calculates the total travel time
    as measured on Earth and onboard the spacecraft.
    """

    # --- Constants and Initial Values ---
    # Using the Fraction type to represent all numbers, as required.
    C_KM_PER_S = Fraction(299792458, 1000)  # Speed of light in km/s
    SECONDS_PER_DAY = Fraction(24 * 60 * 60)
    DAYS_PER_YEAR = Fraction(3652425, 10000)  # Using Gregorian year for accuracy

    # --- Step 1: Calculate Pandora's Recession Velocity ---
    lambda_observed = Fraction(501)
    lambda_emitted = Fraction(500)
    z = (lambda_observed - lambda_emitted) / lambda_emitted
    v_recession = z * C_KM_PER_S
    print(f"Pandora's recession velocity: {float(v_recession):.1f} km/s")


    # --- Step 2: Establish the Mission Target Distance ---
    # This is the key assumption: Initial distance is light-travel time
    # during the spacecraft's 400-day acceleration phase.
    acceleration_duration_days = 400
    D_initial = C_KM_PER_S * acceleration_duration_days * SECONDS_PER_DAY
    print(f"Assumed initial distance to Pandora: {float(D_initial / C_KM_PER_S / SECONDS_PER_DAY):.1f} light-days")

    # --- Utility Function: Sqrt for Fractions ---
    # A custom square root function using Newton's method is required
    # as it's not available on the Wuxing computer.
    def fsqrt(n, iterations=30):
        if n < 0: return Fraction(0)
        if n == 0: return Fraction(0)
        # An initial guess using floating point is efficient.
        x = Fraction(int(math.sqrt(float(n)))) if n > 1 else Fraction(1)
        # Iteratively refine the guess using only basic arithmetic.
        for _ in range(iterations):
            x = (x + n / x) / 2
        return x

    # --- Step 3: Run the Day-by-Day Mission Simulation ---
    v_current = Fraction(40)  # Pioneer's initial velocity in km/s
    distance_closed = Fraction(0)
    onboard_time_days = Fraction(0)
    earth_days_counter = 0

    while distance_closed < D_initial:
        earth_days_counter += 1

        # Update velocity based on the day's acceleration schedule
        if 1 <= earth_days_counter <= 100:
            v_current *= Fraction(104, 100)  # 4.0% relative acceleration
        elif 101 <= earth_days_counter <= 200:
            v_current *= Fraction(102, 100)  # 2.0%
        elif 201 <= earth_days_counter <= 300:
            v_current *= Fraction(101, 100)  # 1.0%
        elif 301 <= earth_days_counter <= 400:
            v_current *= Fraction(1005, 1000) # 0.5%
        # After day 400, velocity remains constant (v_current is not changed)

        # Calculate how much closer Pioneer gets to Pandora today
        distance_closed_today = (v_current - v_recession) * SECONDS_PER_DAY
        distance_closed += distance_closed_today

        # Calculate onboard time using the time dilation formula
        # Δt_onboard = Δt_earth * sqrt(1 - v²/c²)
        v_over_c_squared = (v_current / C_KM_PER_S)**2
        onboard_day_fraction = fsqrt(1 - v_over_c_squared)
        onboard_time_days += onboard_day_fraction

    # --- Step 4: Correct for the final partial day ---
    # The loop overshoots the target; we must correct for this fraction of a day.
    overshoot_dist = distance_closed - D_initial
    last_day_closing_rate = (v_current - v_recession) * SECONDS_PER_DAY
    overshoot_fraction_of_day = overshoot_dist / last_day_closing_rate

    # Correct Earth time
    final_earth_time_days = earth_days_counter - overshoot_fraction_of_day

    # Correct onboard time
    last_day_onboard_frac = fsqrt(1 - (v_current / C_KM_PER_S)**2)
    onboard_time_to_remove = overshoot_fraction_of_day * last_day_onboard_frac
    final_onboard_time_days = onboard_time_days - onboard_time_to_remove

    # --- Final Output ---
    earth_years = final_earth_time_days / DAYS_PER_YEAR
    onboard_years = final_onboard_time_days / DAYS_PER_YEAR
    
    print("\n--- Final Results ---")
    print(f"a (Earth years): {float(earth_years):.1f}")
    print(f"b (Onboard years): {float(onboard_years):.1f}")


solve_pandora_mission()