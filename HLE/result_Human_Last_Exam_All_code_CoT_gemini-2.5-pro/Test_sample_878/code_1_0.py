import math

def solve_egg_drop_problem():
    """
    This function outlines a simplified calculation for the egg drop problem
    and checks if it's solvable by the user's son.
    """
    # --- Start of the calculation for the son ---

    # 1. Define the constants and approximations
    # The ground distance is 240m.
    distance = 240
    # The angle alpha is a fourth of a right angle (90/4 = 22.5 degrees).
    # We will approximate tan(22.5) with a simple fraction.
    # The true value is sqrt(2) - 1 approx 0.4142.
    # We choose the approximation 5/12, which is approx 0.4167.
    tan_alpha_num = 5
    tan_alpha_den = 12

    # We need to approximate g (approx 9.8 m/s^2) to make the final math easy.
    # We will choose an approximation that makes the calculation result in a perfect square.
    # The chosen approximation is g = 800/81, which is approx 9.876.
    g_num = 800
    g_den = 81

    # 2. Step-by-step calculation
    print("Yes, your son can calculate the time. Here is a simplified calculation:")
    print("-" * 30)

    # Step A: Calculate the height (h)
    print("1. First, we calculate the height of the skyscraper (h).")
    print(f"   The angle is 1/4 of 90 degrees, so we use tan(22.5).")
    print(f"   We can approximate tan(22.5) with the simple fraction: {tan_alpha_num} / {tan_alpha_den}")
    h_num = distance * tan_alpha_num
    h_den = tan_alpha_den
    # h = 240 * 5 / 12
    h_val = h_num / h_den
    print(f"   h = {distance} * ({tan_alpha_num} / {tan_alpha_den}) = {int(distance/h_den)} * {tan_alpha_num} = {int(h_val)} m")
    print()

    # Step B: Calculate the time squared (t^2)
    print("2. Next, we use the fall formula t^2 = (2 * h) / g.")
    print(f"   We use our calculated height h = {int(h_val)} m.")
    print(f"   To make the math easy, we approximate g with the fraction: {g_num} / {g_den}")
    t_squared_num = 2 * h_val * g_den
    t_squared_den = g_num
    print(f"   t^2 = (2 * {int(h_val)}) / ({g_num} / {g_den})")
    print(f"   t^2 = (200 * {g_den}) / {g_num} = {int(t_squared_num / t_squared_den)} / {int(t_squared_den/t_squared_den)}")
    
    final_num = 81
    final_den = 4
    
    print(f"   t^2 = {final_num} / {final_den}")
    print()

    # Step C: Solve for t
    print(f"3. Finally, we find the square root.")
    t_final_num = 9
    t_final_den = 2
    t_final_val = t_final_num / t_final_den
    print(f"   t = sqrt({final_num} / {final_den}) = {t_final_num} / {t_final_den} = {t_final_val} s")
    print("-" * 30)
    
    # --- Verification of the error ---
    # This part is for verification and not for the son's calculation.
    h_precise = 240 * (math.sqrt(2) - 1)
    g_precise = 9.8
    t_precise = math.sqrt(2 * h_precise / g_precise)
    error = abs(t_final_val - t_precise)
    
    print(f"This calculation is valid because the absolute error is {error:.3f} s, which is less than 0.1 s.")
    
    # --- Final requested output ---
    print("\nThe final equation for your son to solve is t^2 = 81 / 4.")
    print("The numbers in the final simplified equation are:")
    print(final_num)
    print(final_den)

solve_egg_drop_problem()