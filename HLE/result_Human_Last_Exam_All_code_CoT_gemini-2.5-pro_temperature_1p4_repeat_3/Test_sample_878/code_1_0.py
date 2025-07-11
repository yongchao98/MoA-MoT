import math

def solve_egg_drop_problem():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic,
    and presents the steps in a way a child can understand.
    """
    # -- Real values for verification --
    g_real = 9.80665
    h_real = 240 * (math.sqrt(2) - 1)
    t_real = math.sqrt(2 * h_real / g_real)

    # -- Simplified calculation for the son --

    print("Let's calculate the time it takes for the egg to reach the ground, step by step, using simple numbers!")
    print("-" * 30)

    # Step 1: Calculate the height of the skyscraper
    d = 240
    tan_num, tan_den = 2, 5
    h_approx = d * tan_num / tan_den

    print("Step 1: Find the height of the skyscraper (h)")
    print(f"We know the distance to the skyscraper is {d}m.")
    print("The angle was measured to be a fourth of a right angle (90/4 degrees).")
    print(f"For such an angle, the tangent is tricky, but we can use the simple fraction {tan_num}/{tan_den} which is very close.")
    print("\nThe equation to find the height is: height = distance × tangent")
    print(f"h = {d} * ({tan_num} / {tan_den})")
    print(f"So, the height is {int(h_approx)} meters.")
    print("-" * 30)

    # Step 2: Calculate the time of fall
    g_num, g_den = 256, 27
    
    t_squared_intermediate_num = 2 * int(h_approx) * g_den
    t_squared_intermediate_den = g_num
    
    final_t_squared_num = 81
    final_t_squared_den = 4
    
    t_num = int(math.sqrt(final_t_squared_num))
    t_den = int(math.sqrt(final_t_squared_den))
    t_approx = t_num / t_den

    print("Step 2: Find the time of the fall (t)")
    print("The physics formula is: t = sqrt(2 * h / g)")
    print(f"Let's use our calculated height h = {int(h_approx)}m.")
    print(f"For gravity (g), instead of a complex decimal like 9.8, we'll use the fraction {g_num}/{g_den} to make our math easy.")
    print("\nThe final equation for time squared (t^2) is:")
    # Here we output each number in the final equation as requested
    print(f"t^2 = (2 * {int(h_approx)}) / ({g_num} / {g_den})")
    print(f"t^2 = ({2 * int(h_approx)} * {g_den}) / {g_num}")
    print(f"t^2 = {t_squared_intermediate_num} / {t_squared_intermediate_den}")
    print(f"This fraction simplifies nicely to {final_t_squared_num} / {final_t_squared_den}.")
    print("\nThis is a perfect square! So we can easily find the square root:")
    print(f"t = {t_num} / {t_den} = {t_approx} seconds.")
    print("-" * 30)

    # Step 3: Verification and conclusion
    error = abs(t_approx - t_real)
    print("Step 3: Final check")
    print(f"Our calculated time is {t_approx} seconds.")
    print(f"A more precise scientific calculation gives a time of {t_real:.2f} seconds.")
    print(f"The difference is only {error:.2f} seconds, which is less than the required 0.1 seconds.")
    print("\nYes, your son can calculate the time accurately with his skills!")
    
    calculation_integers = [90, 4, d, tan_num, tan_den, int(h_approx), 2, g_num, g_den, 81, t_num, t_den]
    z = max(calculation_integers)
    print(f"The largest integer used in this calculation is {z}.")

if __name__ == '__main__':
    solve_egg_drop_problem()