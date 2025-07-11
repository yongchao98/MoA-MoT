import math

def solve_egg_drop():
    """
    Calculates the egg drop time using simplified fractions and integers,
    and prints the calculation process.
    """
    # --- Problem constants ---
    distance_d = 240  # meters
    right_angle = 90  # degrees
    angle_divisor = 4

    # --- Approximations for the son ---
    # We need tan(90/4) = tan(22.5). The actual value is sqrt(2) - 1.
    # We approximate sqrt(2) as 17/12 to get a simple fraction.
    # So, tan(22.5) ~= 17/12 - 1 = 5/12.
    tan_alpha_num = 5
    tan_alpha_den = 12

    # We use g = 10 m/s^2 as a simple integer approximation for gravity.
    g_approx = 10

    # --- Calculation with approximations (for the son) ---
    # Height h = d * tan(alpha)
    h_approx = distance_d * tan_alpha_num / tan_alpha_den

    # Time t^2 = 2 * h / g
    t_squared_num = 2 * h_approx
    t_squared_approx = t_squared_num / g_approx

    # --- Error Check (to verify the method) ---
    # More precise calculation
    g_real = 9.8
    h_real = distance_d * (math.sqrt(2) - 1)
    t_real = math.sqrt(2 * h_real / g_real)
    t_approx = math.sqrt(t_squared_approx)
    
    absolute_error = abs(t_real - t_approx)

    if absolute_error < 0.1:
        # The largest integer used in the simplified calculation for the son.
        # The numbers are: 90, 4, 240, 17, 12, 5, 1, 100, 2, 10, 200, 20.
        # The largest is 240.
        answer_code = "Y240"

        print("Yes, your son can calculate the time! Here is a simple way:")
        print("\n--- Step 1: Find the height of the skyscraper (h) ---")
        print(f"First, we find the angle. A right angle is {right_angle} degrees.")
        print(f"The angle alpha is one-fourth of that: {right_angle} / {angle_divisor} degrees.")
        print(f"To find the height, we calculate: h = distance * tan(angle)")
        print(f"A good and simple fraction for tan({right_angle}/{angle_divisor}) is {tan_alpha_num}/{tan_alpha_den}.")
        print(f"So, h = {distance_d} * {tan_alpha_num} / {tan_alpha_den} = {int(h_approx)} meters.")

        print("\n--- Step 2: Find the time (t) ---")
        print(f"The physics formula is: t * t = 2 * h / g")
        print(f"We can use a simple value for gravity, g = {g_approx}.")
        print(f"So, t * t = (2 * {int(h_approx)}) / {g_approx}")
        print(f"t * t = {int(t_squared_num)} / {g_approx}")
        print(f"t * t = {int(t_squared_approx)}")
        print(f"\nThe time is the square root of {int(t_squared_approx)}, which is about {t_approx:.1f} seconds.")
    else:
        answer_code = "N0"
        print("No, the approximations lead to an error of more than 0.1 seconds.")

# Execute the function
solve_egg_drop()