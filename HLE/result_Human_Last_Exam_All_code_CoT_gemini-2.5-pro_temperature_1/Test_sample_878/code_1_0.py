import math

def solve_egg_drop_problem():
    """
    Calculates the egg drop time using simplified fractions and checks the accuracy.
    It prints the calculation steps in a way a child can follow.
    """
    # --- Part 1: Define the problem with precise numbers ---
    distance = 240  # meters
    g_actual = 9.8  # m/s^2
    # height = distance * tan(22.5 degrees) = distance * (sqrt(2) - 1)
    height_actual = distance * (math.sqrt(2) - 1)
    time_actual = math.sqrt(2 * height_actual / g_actual)

    # --- Part 2: Find "clever" fractional approximations ---
    # We want to approximate g and sqrt(2) so that the final calculation is simple.
    # We choose g = 49/5, which is exactly 9.8.
    g_num, g_den = 49, 5
    # We need to approximate sqrt(2). Let's try sqrt(2) ~= 415/294.
    # This choice is made specifically to make the final square root easy.
    sqrt2_num, sqrt2_den = 415, 294

    # --- Part 3: Perform the calculation as the son would ---
    # The numbers used to set up the calculation are:
    # 2 (from the formula), distance (240), sqrt2_num (415), sqrt2_den (294),
    # g_num (49), g_den (5)
    
    # Let's show the steps:
    print("Yes, your son can solve this! Here is a simple way to do the calculation.")
    print("-" * 60)
    print("The goal is to calculate the time 't'. The formula is:")
    print("Time = sqrt( (2 * Height) / g )")
    print("\nFirst, we find the Height of the skyscraper:")
    print("Height = Distance * (sqrt(2) - 1)")
    print(f"The distance is {distance} m.")
    print(f"We will use the fraction {sqrt2_num}/{sqrt2_den} as an approximation for sqrt(2).")
    print(f"So, sqrt(2) - 1 is approximately {sqrt2_num}/{sqrt2_den} - 1 = {sqrt2_num - sqrt2_den}/{sqrt2_den} = 121/294.")
    
    print("\nNext, we need a value for 'g', the acceleration of gravity.")
    print(f"We will use the fraction {g_num}/{g_den}, which is exactly 9.8.")
    
    print("\nNow let's put all the numbers into the time formula:")
    print(f"Time = sqrt( (2 * {distance} * (121 / {sqrt2_den})) / ({g_num} / {g_den}) )")
    print("\nThis looks complicated, but let's simplify it. First, multiply all the numbers on top and all the numbers on the bottom:")
    print(f"Time = sqrt( (2 * {distance} * 121 * {g_den}) / ({sqrt2_den} * {g_num}) )")
    
    # Perform the simplification
    # t^2 = (2 * 240 * 121 * 5) / (294 * 49)
    # t^2 = (2400 * 121) / 14406
    # We simplify by noticing 294 = 6 * 49, and 240 * 5 / (6) = 200 * 2 / 2 = 400/3 *
    # Let's do it differently: 294 = 2*3*49
    # t^2 = (2 * 240 * 121 * 5) / (2 * 3 * 49 * 49)
    # t^2 = (240 * 121 * 5) / (3 * 49 * 49)
    # t^2 = (80 * 121 * 5) / (49 * 49)
    # t^2 = (400 * 121) / (49 * 49)
    # t^2 = (20*20 * 11*11) / (49*49)
    
    print("\nAfter some simplification, the fraction becomes much nicer:")
    print("Time = sqrt( (400 * 121) / (49 * 49) )")
    print("\nThis is the same as:")
    print("Time = sqrt( (20 * 20 * 11 * 11) / (49 * 49) )")
    print("\nNow the square root is easy!")
    print("Time = (20 * 11) / 49")
    final_num = 20 * 11
    final_den = 49
    print(f"Time = {final_num} / {final_den}")

    time_approx = final_num / final_den
    error = abs(time_approx - time_actual)

    print("-" * 60)
    print(f"\nThe final answer is {time_approx:.2f} seconds.")
    print(f"The more precise answer is {time_actual:.2f} seconds.")
    print(f"The error is only {error:.3f} seconds, which is less than 0.1 seconds. Great success!")

solve_egg_drop_problem()

# The largest integer used in the setup of the calculation was 415.
# The calculation is possible with the required accuracy.
# So the final answer is Y415.
print("\n<<<Y415>>>")