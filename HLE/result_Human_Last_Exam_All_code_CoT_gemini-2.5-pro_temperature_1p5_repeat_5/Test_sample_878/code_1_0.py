import math

def solve_egg_drop_problem():
    """
    This function calculates the egg drop time using approximations
    with small integers and fractions, then prints the calculation steps.
    """

    # Step 1: Define the constants and approximations
    # Distance from the skyscraper
    d = 240  # meters
    # Angle alpha is a fourth of a right angle (90 degrees), so alpha = 22.5 degrees.
    # The tangent of 22.5 degrees is sqrt(2) - 1.
    # We can approximate sqrt(2) as 7/5.
    # So, tan(22.5) ≈ 7/5 - 1 = 2/5.
    tan_alpha_num = 2
    tan_alpha_den = 5
    
    # We use a fractional approximation for gravity, g.
    # g ≈ 9.8 m/s^2 = 98/10 = 49/5 m/s^2.
    g_num = 49
    g_den = 5

    # Step 2: Calculate the height 'h' using our approximation for tan(alpha)
    h_val = d * tan_alpha_num / tan_alpha_den

    # Step 3: Calculate t^2 = 2 * h / g
    t_squared_num = 2 * h_val * g_den
    t_squared_den = g_num
    
    # Step 4: The numerator is 960, which is not a perfect square.
    # Let's find the nearest perfect square to make it easy for a child.
    # 30*30 = 900 and 31*31 = 961. 960 is very close to 961.
    t_squared_num_approx = 961
    final_t_num = int(math.sqrt(t_squared_num_approx))
    final_t_den = int(math.sqrt(t_squared_den))
    
    # Let's present the logic and the calculation.
    print("Yes, your son can calculate the time! Here is a simple way:")
    print("\nThe formula to calculate the fall time 't' from height 'h' is: t^2 = (2 * h) / g")
    print("\nFirst, we find the height 'h'.")
    print(f"h = distance * tan(angle) = {d} * tan(22.5°)")
    print(f"We can use a simple fraction for tan(22.5°), which is {tan_alpha_num}/{tan_alpha_den}.")
    print(f"So, h = {d} * {tan_alpha_num} / {tan_alpha_den} = {int(h_val)} meters.")
    
    print("\nNext, we use a good fraction for gravity, g.")
    print(f"g ≈ {g_num}/{g_den} m/s².")
    
    print("\nNow, we put all the numbers into the time formula:")
    print(f"t^2 = (2 * {int(h_val)}) / ({g_num}/{g_den})")
    print(f"t^2 = (2 * {int(h_val)} * {g_den}) / {g_num} = {int(t_squared_num)} / {t_squared_den}")

    print("\nTo find 't', we need to take the square root. Let's make it easier!")
    print(f"The number {int(t_squared_num)} is very close to {t_squared_num_approx}, which is a perfect square: {final_t_num} * {final_t_num}.")
    print(f"So we can say: t^2 ≈ {t_squared_num_approx} / {t_squared_den}")
    
    print("\nTaking the square root of the top and bottom gives us the final answer:")
    print(f"t ≈ sqrt({t_squared_num_approx}) / sqrt({t_squared_den}) = {final_t_num} / {final_t_den} seconds.")

solve_egg_drop_problem()