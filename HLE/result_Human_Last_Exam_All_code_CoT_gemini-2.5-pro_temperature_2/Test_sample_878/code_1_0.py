import math

def solve_egg_drop():
    """
    Calculates the time for an egg to fall from a skyscraper using simple integer
    and fraction arithmetic, suitable for a young learner.
    """

    print("Yes, the calculation is possible with an error less than 0.1s.")
    print("Here is a simple calculation using small integers and fractions:")
    print("-----------------------------------------------------------------")

    # Step 1: Find the height (h) of the skyscraper.
    # The viewing angle is a fourth of a right angle (90 degrees).
    angle_total = 90
    angle_divisor = 4
    angle = angle_total / angle_divisor

    # The distance from the building is 240m.
    distance = 240
    
    # We approximate tan(22.5) ≈ 5/12 to keep numbers simple.
    # tan(22.5) is exactly sqrt(2)-1 ≈ 0.414, and 5/12 ≈ 0.417, which is close enough.
    tan_num = 5
    tan_den = 12
    
    print(f"Step 1: First, let's find the skyscraper's height (h).")
    print(f"The angle is {angle_total} / {angle_divisor} = {angle} degrees.")
    print(f"The distance to the skyscraper is {distance}m.")
    print(f"We can find the height using: h = {distance} * tan({angle})")
    print(f"To make it simple, we can use the fraction {tan_num}/{tan_den} for tan({angle}).")
    print(f"The equation for the height is: h = {distance} * {tan_num} / {tan_den}")

    height = int(distance * tan_num / tan_den)
    print(f"So, the height is {height}m.")
    print("-----------------------------------------------------------------")
    
    # Step 2: Find the time (t) it takes for the egg to fall.
    # The formula is t = sqrt(2 * h / g).
    # We use g ≈ 10 m/s^2 for simplicity (actual value is ~9.8).
    g_approx = 10
    factor_2 = 2

    print(f"Step 2: Now, let's calculate the fall time (t).")
    print(f"The physics equation is: t = sqrt( (2 * h) / g )")
    print(f"We use our height h = {height}m and g = {g_approx} m/s^2 to make the math easy.")
    print(f"The equation for time is: t = sqrt( ({factor_2} * {height}) / {g_approx} )")
    
    time_squared = int((factor_2 * height) / g_approx)
    print(f"t = sqrt({time_squared})")
    print("-----------------------------------------------------------------")

    # Step 3: Approximate the square root.
    print(f"Step 3: Finally, let's estimate sqrt({time_squared}).")
    factor_a = 4
    factor_b = 5
    print(f"We can rewrite {time_squared} as {factor_a} * {factor_b}, so sqrt({time_squared}) = sqrt({factor_a}) * sqrt({factor_b}) = 2 * sqrt({factor_b}).")
    
    # We approximate sqrt(5) ≈ 9/4, because (9/4)^2 = 81/16 = 5.0625, which is very close to 5.
    sqrt5_num = 9
    sqrt5_den = 4
    print(f"A good and simple fraction for sqrt({factor_b}) is {sqrt5_num}/{sqrt5_den}.")
    print(f"Our final equation is: t = 2 * {sqrt5_num} / {sqrt5_den}")
    
    final_time_num = factor_2 * sqrt5_num
    final_time_den = sqrt5_den
    final_time = final_time_num / final_time_den
    
    print(f"The final answer is: t = {final_time_num} / {final_time_den} = {final_time} seconds.")
    
    # The question asks to answer Yz with z being the largest integer in the calculation.
    # The integers used in the printed equations are: 90, 4, 22.5 (not int), 240, 5, 12,
    # 100, 2, 10, 20. The largest integer is 240.
    largest_integer_used = 240
    print(f"<<<Y{largest_integer_used}>>>")

solve_egg_drop()