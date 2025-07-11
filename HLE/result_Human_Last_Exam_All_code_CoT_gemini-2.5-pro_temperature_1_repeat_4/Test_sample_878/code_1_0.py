import math

def solve_egg_drop_problem():
    """
    Calculates the time for an egg to fall from a skyscraper using simplified
    fractions for a child, and determines if the result is accurate enough.
    """
    # Given values from the problem
    distance = 240  # meters

    # Approximations chosen to simplify the calculation for a child
    # tan(22.5°) is approximated as 27/64
    tan_alpha_num = 27
    tan_alpha_den = 64
    # g is approximated as 10 m/s^2
    g_approx = 10

    # --- Calculation Explained for the Son ---
    print("Let's calculate the time it takes for the egg to reach the ground!")
    print("\n--- Step 1: How tall is the skyscraper? ---")
    print(f"We are standing {distance} meters away from the skyscraper.")
    print("The angle to the top is a quarter of a right angle (90 / 4 = 22.5 degrees).")
    print("To find the height, we need a magic fraction for this angle. Let's use a close one: 27 / 64.")
    print(f"Height = {distance} * ( {tan_alpha_num} / {tan_alpha_den} )")
    
    # Calculate height h = 240 * 27 / 64 = 405 / 4
    h_num = distance * tan_alpha_num
    h_den = tan_alpha_den
    print(f"Height = {h_num} / {h_den} meters. Let's simplify this to 405 / 4 meters.")

    print("\n--- Step 2: How long does the egg fall? ---")
    print("We have a formula for falling time: Time² = (2 * Height) / g")
    print("For gravity 'g', let's use the simple number 10.")

    print("\nNow we put all our numbers into the final equation:")
    print("Time² = (2 * (405 / 4)) / 10")
    
    # Simplify the equation step-by-step
    print("Time² = (810 / 4) / 10")
    print("Time² = 810 / 40")
    print("Time² = 81 / 4")
    
    print("\nTo find the time, we need a number that, when multiplied by itself, gives 81 / 4.")
    print("Since 9 * 9 = 81, and 2 * 2 = 4, our answer is:")
    print("Time = 9 / 2 = 4.5 seconds.")

    # --- Verification of Accuracy (for the parent) ---
    h_actual = 240 * (math.sqrt(2) - 1)
    t_actual = math.sqrt(2 * h_actual / 9.8)
    error = abs(4.5 - t_actual)

    print("\n--- Checking Our Work ---")
    print(f"The actual time is about {t_actual:.3f} seconds.")
    print(f"Our calculated time is 4.5 seconds.")
    print(f"The difference (error) is only {error:.3f} seconds, which is less than 0.1 seconds.")
    print("So, YES, your son can calculate the time with great accuracy using these simple numbers!")

    # --- Finding the largest integer 'z' ---
    # The initial calculation is Time^2 = (2 * distance * tan_approx) / g_approx
    # The integers involved are 2, distance (240), tan_num (27), tan_den (64), g (10)
    largest_integer = max(2, distance, tan_alpha_num, tan_alpha_den, g_approx)
    print(f"\nThe integers we used to set up the final equation were 2, {distance}, {tan_alpha_num}, {tan_alpha_den}, and {g_approx}.")
    print(f"The largest integer in this calculation is {largest_integer}.")

solve_egg_drop_problem()
