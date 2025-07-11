import math

def solve_egg_drop():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic,
    and prints the steps in a way a child can understand.
    """

    print("Hello! Let's figure out how long it takes for the egg to drop. It's a fun two-step problem!")
    print("\n--- Step 1: Find the height of the skyscraper ---")
    
    # h = 240 * tan(22.5)
    # tan(22.5) = sqrt(2) - 1 ≈ 0.414. A good fraction is 5/12 = 0.4166.
    # g ≈ 10 m/s^2
    # sqrt(5) ≈ 9/4 = 2.25
    
    distance = 240
    tan_alpha_num, tan_alpha_den = 5, 12
    
    print(f"We are {distance} meters away from the skyscraper.")
    print("The angle to the top is special. To find the height (h), we need to multiply our distance by the tangent of that angle.")
    print(f"A good and simple fraction for the tangent of our angle is {tan_alpha_num}/{tan_alpha_den}.")
    print(f"\nSo, the height calculation is:")
    print(f"h = {distance} * ({tan_alpha_num} / {tan_alpha_den})")
    
    h_intermediate = distance // tan_alpha_den
    height = h_intermediate * tan_alpha_num
    
    print(f"h = {h_intermediate} * {tan_alpha_num}")
    print(f"h = {height} meters")
    
    print("\n--- Step 2: Calculate the falling time ---")
    
    g_approx = 10
    
    print(f"The skyscraper is {height} meters tall. Wow!")
    print(f"The formula for the fall time (t) is t = square root of (2 * h / g).")
    print(f"For gravity (g), we can use the simple number {g_approx} instead of 9.8.")
    
    t_squared_num = 2 * height
    t_squared = t_squared_num / g_approx
    
    print("\nSo, the time squared calculation is:")
    print(f"t*t = (2 * {height}) / {g_approx}")
    print(f"t*t = {t_squared_num} / {g_approx}")
    print(f"t*t = {int(t_squared)}")
    
    print("\nNow we need to find the square root of 20.")
    print("t = square root of (4 * 5) = 2 * (square root of 5)")
    print("A good fraction for the square root of 5 is 9/4.")

    # Final calculation t ≈ 2 * (9/4)
    final_calc_factor1 = 2
    final_calc_num, final_calc_den = 9, 4
    final_time = (final_calc_factor1 * final_calc_num) / final_calc_den

    print("\n--- The Final Answer ---")
    print("Our final equation to get the time in seconds is:")
    print(f"t = {final_calc_factor1} * {final_calc_num} / {final_calc_den}")
    print("\nThe numbers in this final equation are:")
    print(final_calc_factor1)
    print(final_calc_num)
    print(final_calc_den)
    
    print(f"\nSo, the time is {final_calc_factor1*final_calc_num}/{final_calc_den} seconds, which is {final_time} seconds!")
    
    # Verification
    # True value: t = sqrt(2 * 240 * (sqrt(2)-1) / 9.8) ≈ 4.504s
    # Our value: 4.5s
    # Absolute error: |4.504 - 4.5| = 0.004s, which is less than 0.1s.
    # Largest integer in calculation: distance = 240.
    
    print("\nThis calculation is accurate enough for your son!")
    print("<<<Y240>>>")

solve_egg_drop()