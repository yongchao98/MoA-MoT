import math

# Wuxing C Program Simulation

# This script simulates the calculations that would be performed by a C program
# on the hypothetical Wuxing computer to solve the Pandora probe problem.
# The Wuxing computer uses a special 'frac' type for floating-point arithmetic.
# This simulation uses standard Python floats for simplicity, but follows the
# analytical steps that would be necessary on the constrained hardware.

def solve_pandora_mission():
    """
    Calculates the travel time to Pandora for both Earth and Pioneer reference frames.
    """

    # --- Step 1: Define constants and physical parameters ---

    # Speed of light in km/s
    c = 300000.0

    # Pioneer's initial speed in km/s
    v_initial = 40.0

    # Pandora's recessional velocity calculation
    # Redshift z = (lambda_obs - lambda_emit) / lambda_emit
    z = (501.0 - 500.0) / 500.0
    # Recessional velocity v = z * c (for non-relativistic speeds)
    v_pandora = z * c

    # Distance to Pandora calculation using Hubble's Law (v = H0 * D)
    # We assume the standard Hubble Constant H0 ≈ 70 km/s/Mpc
    H0 = 70.0  # km/s per Megaparsec
    # Conversion factor for Megaparsec to km
    Mpc_to_km = 3.0857e19
    # Initial distance to Pandora D = v_pandora / H0
    D_initial = (v_pandora / H0) * Mpc_to_km

    # Time conversion constants
    secs_per_day = 86400.0
    days_per_year = 365.25

    # --- Step 2: Simulate the 400-day acceleration phase ---

    # This part calculates the probe's final velocity after 400 days
    # and the total time elapsed on the probe's clock during this phase.

    v = v_initial
    total_onboard_secs_accel = 0.0

    # A helper function to calculate square roots, as Wuxing has no built-in sqrt.
    # This uses the Babylonian method, which only requires basic arithmetic.
    def numerical_sqrt(s, iterations=20):
        if s < 0: return 0
        if s == 0: return 0
        # Initial guess for the square root
        x = 1.0
        for _ in range(iterations):
            x = 0.5 * (x + s / x)
        return x

    for day in range(400):
        # Calculate time dilation for the current day
        # Lorentz factor gamma = 1 / sqrt(1 - v^2/c^2)
        # Onboard time delta_tau = delta_t / gamma = delta_t * sqrt(1 - v^2/c^2)
        if v >= c:
            v = c
            gamma_inv = 0
        else:
            gamma_inv = numerical_sqrt(1.0 - (v/c)**2)
        
        delta_tau_secs = secs_per_day * gamma_inv
        total_onboard_secs_accel += delta_tau_secs

        # Update velocity for the next day based on the specified acceleration schedule
        if day < 100:
            v *= 1.04  # 4% acceleration
        elif day < 200:
            v *= 1.02  # Acceleration halved to 2%
        elif day < 300:
            v *= 1.01  # Acceleration halved to 1%
        elif day < 400:
            v *= 1.005 # Acceleration halved to 0.5%
        # After day 399, acceleration stops.

    v_final = v

    # --- Step 3: Calculate the constant velocity cruise phase ---

    # The distance traveled during the 400-day acceleration is negligible
    # compared to the total interstellar distance, so we can simplify the calculation.
    # We calculate the time to cover the initial distance D_initial at a
    # constant closing speed.

    # Closing speed = Pioneer's final speed - Pandora's recessional speed
    v_closing = v_final - v_pandora

    # Time for the cruise phase, as measured on Earth
    T_cruise_secs_earth = D_initial / v_closing

    # --- Step 4: Calculate total journey times in years ---

    # a. Total time as measured on Earth
    # Total time = acceleration phase time + cruise phase time
    T_total_secs_earth = (400.0 * secs_per_day) + T_cruise_secs_earth
    T_total_years_earth = T_total_secs_earth / (secs_per_day * days_per_year)

    # b. Total time as measured by Pioneer's onboard clock
    # First, find the onboard time for the cruise phase by applying time dilation
    gamma_final_inv = numerical_sqrt(1.0 - (v_final/c)**2)
    T_cruise_secs_onboard = T_cruise_secs_earth * gamma_final_inv

    # Total onboard time = onboard time from accel phase + onboard time from cruise phase
    T_total_secs_onboard = total_onboard_secs_accel + T_cruise_secs_onboard
    T_total_years_onboard = T_total_secs_onboard / (secs_per_day * days_per_year)

    # --- Step 5: Output the final answer ---
    # The problem asks for the answer in the format a:b, where a and b are the
    # integer results. We round the final values to the nearest integer.

    answer_a = int(round(T_total_years_earth))
    answer_b = int(round(T_total_years_onboard))

    # The prompt also says: "Remember in the final code you still need to output each number in the final equation!"
    # This is interpreted as printing the final two numbers that form the answer.
    print(f"{answer_a}:{answer_b}")

solve_pandora_mission()