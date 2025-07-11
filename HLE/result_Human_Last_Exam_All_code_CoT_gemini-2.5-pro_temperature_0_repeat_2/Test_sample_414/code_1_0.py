import math

def solve_ingrowth_time():
    """
    Calculates the time between chemical separation and the first measurement
    based on the ingrowth of Yttrium-90 from Strontium-90.
    """
    # --- Given Data ---
    # Initial measured activity in kBq/mL
    activity_1 = 1.4
    # Activity measured 14 days later in kBq/mL
    activity_2 = 2.1
    # Time between the two measurements in days
    delta_t_days = 14.0
    # Half-life of Yttrium-90 in hours
    t_half_y90_hours = 64.1

    # --- Calculations ---
    # Convert half-life to days
    t_half_y90_days = t_half_y90_hours / 2.4e1

    # Calculate the decay constant (lambda) for Y-90 in days^-1
    # lambda = ln(2) / T_half
    lambda_y = math.log(2) / t_half_y90_days

    # The ingrowth of Y-90 activity (A_y) follows the equation:
    # A_y(t) = A_sr * (1 - exp(-lambda_y * t))
    # We have two points in time, t1 and t1 + delta_t.
    # A1 = A_sr * (1 - exp(-lambda_y * t1))
    # A2 = A_sr * (1 - exp(-lambda_y * (t1 + delta_t)))
    # By dividing the two equations, we eliminate the unknown A_sr:
    # A2 / A1 = (1 - exp(-lambda_y * (t1 + delta_t))) / (1 - exp(-lambda_y * t1))

    # Let's solve for t1. Let x = exp(-lambda_y * t1)
    # ratio = (1 - x * exp(-lambda_y * delta_t)) / (1 - x)
    # After rearranging, we get:
    # x = (ratio - 1) / (ratio - exp(-lambda_y * delta_t))
    # And t1 = -ln(x) / lambda_y

    # Calculate the ratio of the activities
    ratio = activity_2 / activity_1

    # Calculate the decay factor for the 14-day period
    decay_factor = math.exp(-lambda_y * delta_t_days)

    # Calculate the intermediate term 'x'
    x = (ratio - 1) / (ratio - decay_factor)

    # Calculate t1, the time from separation to the first measurement
    t1_days = -math.log(x) / lambda_y

    # --- Output Results ---
    print("This script calculates the time elapsed between a chemical separation and a radioactivity measurement.")
    print("The calculation is based on the ingrowth of a daughter nuclide (Y-90) from a parent (Sr-90).\n")
    print("--- Key Parameters ---")
    print(f"Half-life of Y-90: {t_half_y90_hours} hours")
    print(f"Time between measurements: {delta_t_days} days")
    print(f"First activity (A1): {activity_1} kBq/mL")
    print(f"Second activity (A2): {activity_2} kBq/mL\n")

    print("--- Calculation Steps ---")
    print("The final time `t1` is calculated using the formula: t1 = -ln(x) / lambda_y")
    print("The numbers used in this final equation are:")
    print(f"  Decay constant (lambda_y): {lambda_y:.6f} days^-1")
    print(f"  Intermediate variable (x): {x:.6f}")
    print("  (where x is derived from the activity ratio and decay over 14 days)\n")

    print("--- Final Answer ---")
    print("The approximate time between the chemical separation and the first analysis is:")
    print(f"{t1_days:.2f} days")

solve_ingrowth_time()
<<<4.16>>>