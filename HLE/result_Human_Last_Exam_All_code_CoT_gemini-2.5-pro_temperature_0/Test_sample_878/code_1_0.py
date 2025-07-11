import math

def solve_and_explain_egg_drop():
    """
    This function provides a simplified step-by-step calculation for the egg drop problem,
    suitable for a young learner. It uses fractional approximations to keep numbers simple
    and ensures the final error is within the specified tolerance.
    """

    # --- Define problem parameters and approximations ---
    distance = 240  # meters

    # Approximation for tan(22.5 degrees), which is ~0.4142
    tan_approx_num, tan_approx_den = 5, 12

    # Approximation for g, which is ~9.8 m/s^2
    g_approx_num, g_approx_den = 49, 5

    # --- Step 1: Calculate the approximate height (h) ---
    # h = 240 * (5/12) = 100
    height_approx = (distance * tan_approx_num) // tan_approx_den

    # --- Step 2: Calculate the approximate time squared (t^2) ---
    # t^2 = (2 * h) / g = (2 * 100) / (49/5) = 1000 / 49
    t_squared_num = 2 * height_approx * g_approx_den
    t_squared_den = g_approx_num

    # --- Step 3: Calculate the approximate time (t) ---
    # t = sqrt(1000/49) = 10 * sqrt(10) / 7
    # We need to approximate sqrt(10), which is ~3.162
    # A good simple fraction is 16/5 = 3.2
    sqrt10_approx_num, sqrt10_approx_den = 16, 5
    # t = (10 * (16/5)) / 7 = 32 / 7
    final_t_num = (10 * sqrt10_approx_num) // sqrt10_approx_den
    final_t_den = 7

    # --- Verification (for our reference, not shown to the son) ---
    t_approx_val = final_t_num / final_t_den
    t_true_val = math.sqrt(2 * distance * math.tan(math.radians(22.5)) / 9.80665)
    is_accurate_enough = abs(t_approx_val - t_true_val) < 0.1

    if not is_accurate_enough:
        print("N0")
        return

    # --- Present the calculation ---
    print("Yes, your son can calculate the time! Here is a simple way to do it:\n")

    print("Step 1: Find the height of the skyscraper (h).")
    print(f"The angle is one-fourth of a right angle, which is 22.5 degrees.")
    print(f"We can use the simple fraction {tan_approx_num}/{tan_approx_den} for tan(22.5).")
    print(f"h = {distance} * {tan_approx_num} / {tan_approx_den} = {height_approx} meters\n")

    print("Step 2: Set up the equation for the fall time (t).")
    print("The formula is t^2 = 2 * h / g.")
    print(f"We can use g = {g_approx_num}/{g_approx_den}, which is 9.8.")
    print(f"t^2 = (2 * {height_approx} * {g_approx_den}) / {g_approx_num}")
    print(f"t^2 = {t_squared_num} / {t_squared_den}\n")

    print("Step 3: Find the final answer for t.")
    print(f"To find t, we take the square root: t = sqrt({t_squared_num}) / sqrt({t_squared_den}) = (10 * sqrt(10)) / 7.")
    print(f"A good fraction for sqrt(10) is {sqrt10_approx_num}/{sqrt10_approx_den}.")
    print(f"So, t = (10 * {sqrt10_approx_num}/{sqrt10_approx_den}) / {final_t_den}")
    
    final_equation = f"t = {final_t_num} / {final_t_den}"
    print(f"The final equation is: {final_equation}")
    
    # As requested: "output each number in the final equation"
    print(f"The numbers in the final equation are {final_t_num} and {final_t_den}.")

solve_and_explain_egg_drop()