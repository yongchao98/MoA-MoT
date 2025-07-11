import math
from fractions import Fraction

def solve_pandora_mission():
    """
    Simulates the Pioneer mission to determine the travel time in Earth years
    and onboard years, as calculated by the hypothetical Wuxing computer.
    """
    # --- Part 1: Constants and Initial Values ---
    # Using Fraction to emulate the Wuxing computer's `frac` type
    C = Fraction(300000)          # Speed of light in km/s
    V0 = Fraction(40)             # Pioneer's initial velocity in km/s
    SECONDS_PER_DAY = Fraction(86400)
    YEAR_DIVISOR = 365            # Integer years as would be used in a simple C program

    # Daily velocity multipliers
    MULT1 = Fraction(104, 100)
    MULT2 = Fraction(102, 100)
    MULT3 = Fraction(101, 100)
    MULT4 = Fraction(1005, 1000)

    # --- Part 2: Calculate Pandora's Recession Velocity ---
    z = Fraction(501 - 500, 500) # Redshift z = 1/500
    v_pandora = z * C            # v = z*c, so 600 km/s
    print(f"Pandora's recession velocity calculated as {v_pandora.numerator}/{v_pandora.denominator} km/s.")

    # --- Part 3: Simulation Loop ---
    day = 0
    current_velocity = V0
    total_distance = Fraction(0)
    total_time_correction_sec = Fraction(0) # Time lost due to dilation

    while True:
        day += 1

        # Distance traveled on this specific day
        distance_this_day = current_velocity * SECONDS_PER_DAY
        total_distance += distance_this_day
        
        # Calculate time dilation for this day using the approximation: dt' = dt * (1 - 0.5*(v/c)^2)
        beta_sq = (current_velocity / C)**2
        correction_this_day = Fraction(1, 2) * beta_sq * SECONDS_PER_DAY
        total_time_correction_sec += correction_this_day

        # Earth time elapsed so far
        total_earth_time_sec = day * SECONDS_PER_DAY
        
        # Check for rendezvous: Pioneer's avg speed must equal Pandora's speed.
        # total_distance / total_earth_time_sec >= v_pandora
        if total_distance >= v_pandora * total_earth_time_sec:
            break
            
        # Update Pioneer's velocity for the *next* day
        if day < 100:
            current_velocity *= MULT1
        elif day < 200:
            current_velocity *= MULT2
        elif day < 300:
            current_velocity *= MULT3
        elif day < 400:
            current_velocity *= MULT4
        # After 400 days, acceleration stops and velocity remains constant.
    
    # --- Part 4: Calculate Final Answers in Years ---
    earth_days_to_arrival = day
    onboard_time_sec = (earth_days_to_arrival * SECONDS_PER_DAY) - total_time_correction_sec
    onboard_days_to_arrival = onboard_time_sec / SECONDS_PER_DAY
    
    # A Wuxing C program would use integer division for years
    earth_years_answer = earth_days_to_arrival // YEAR_DIVISOR
    # Onboard days must be converted to an integer before division, simulating a cast
    onboard_years_answer = math.floor(onboard_days_to_arrival) // YEAR_DIVISOR

    print("\n--- Final Calculation ---")
    print(f"a. Earth time to arrival is {earth_days_to_arrival} days.")
    print(f"Final equation for Earth years (a): {earth_days_to_arrival} / {YEAR_DIVISOR} = {earth_years_answer}")
    
    print(f"\nb. Onboard time to arrival is {float(onboard_days_to_arrival):.4f} days.")
    print(f"Final equation for onboard years (b): floor({float(onboard_days_to_arrival):.4f}) / {YEAR_DIVISOR} = {math.floor(onboard_days_to_arrival)} / {YEAR_DIVISOR} = {onboard_years_answer}")

    print("\n--- Program Output ---")
    print(f"{earth_years_answer}:{onboard_years_answer}")

solve_pandora_mission()