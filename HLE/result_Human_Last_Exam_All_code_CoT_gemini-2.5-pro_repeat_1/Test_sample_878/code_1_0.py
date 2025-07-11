import math

def solve_egg_drop():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic.
    """
    # Step 1: Define initial values and approximations
    dist_from_building = 240
    g_approx = 10
    tan_alpha_num = 27
    tan_alpha_den = 64

    print("Yes, the calculation is possible. Here is a simplified method for your son:")
    print("-" * 50)
    print("Step 1: Calculate the height 'h' of the skyscraper.")
    print(f"The height is the distance ({dist_from_building}m) times tan(22.5 degrees).")
    print(f"We approximate tan(22.5 degrees) as the fraction {tan_alpha_num}/{tan_alpha_den}.")
    print(f"h = {dist_from_building} * ({tan_alpha_num} / {tan_alpha_den})")
    
    # h = 240 * 27 / 64. Simplify by dividing 240 and 64 by 8.
    # h = (30 * 27) / 8. Simplify by dividing 30 and 8 by 2.
    # h = (15 * 27) / 4 = 405 / 4
    h_simplified_num = 405
    h_simplified_den = 4
    print(f"h = {h_simplified_num} / {h_simplified_den} meters.")
    print("-" * 50)

    print("Step 2: Calculate the time 't' to fall.")
    print("The formula is t^2 = (2 * h) / g.")
    print(f"We use our calculated h and approximate g as {g_approx}.")
    print(f"t^2 = (2 * ({h_simplified_num}/{h_simplified_den})) / {g_approx}")
    
    # t^2 = (2 * 405/4) / 10 = (405/2) / 10 = 405 / 20
    # Simplify 405/20 by dividing both by 5
    t_sq_final_num = 81
    t_sq_final_den = 4
    print(f"t^2 = {h_simplified_num * 2} / {h_simplified_den * g_approx} = {h_simplified_num} / {(h_simplified_den//2) * g_approx} = {h_simplified_num} / {2*g_approx} = 405 / 20")
    print(f"Simplifying this fraction gives: t^2 = {t_sq_final_num} / {t_sq_final_den}")
    print("-" * 50)
    
    # Another way to show the calculation leading to the largest integer:
    # t^2 = (2 * 240 * 27/64) / 10 = (48 * 27) / 64 = 1296 / 64
    # The largest integer appearing in this calculation path is 1296.
    
    print("Step 3: Find the final time 't' by taking the square root.")
    t_final_num = int(math.sqrt(t_sq_final_num))
    t_final_den = int(math.sqrt(t_sq_final_den))
    print(f"t = sqrt({t_sq_final_num}) / sqrt({t_sq_final_den})")
    print(f"The final equation is: t = {t_final_num} / {t_final_den}")

    # Accuracy check
    t_approx = t_final_num / t_final_den
    t_real = math.sqrt(2 * 240 * (math.sqrt(2) - 1) / 9.8)
    error = abs(t_real - t_approx)
    print(f"\nOur calculated time is {t_approx} seconds.")
    print(f"The actual time is about {t_real:.4f} seconds.")
    print(f"The error of {error:.4f} seconds is less than 0.1 seconds, so this method works!")

solve_egg_drop()