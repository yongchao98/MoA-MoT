import math

def calculate_and_explain():
    """
    Calculates the egg drop time using simple integer and fraction arithmetic,
    and prints the step-by-step process.
    """

    # --- Values from the problem ---
    distance = 240

    # --- Approximations for the calculation ---
    # We will approximate sqrt(2) as 17/12 to calculate tan(22.5)
    sqrt2_num, sqrt2_den = 17, 12
    # We will approximate g as 10 m/s^2
    g = 10
    # We will approximate sqrt(5) as 9/4
    sqrt5_num, sqrt5_den = 9, 4

    # List to keep track of all integers used in the explanation
    integers_used = []

    # --- Part 1: Calculate the height (h) ---
    print("Yes, your son can calculate the time. Here is a possible calculation:")
    print("-" * 30)
    print("Part 1: Calculate the height of the skyscraper (h)")
    print(f"The distance 'd' from the skyscraper is {distance}m.")
    integers_used.append(distance)

    # tan(22.5) = sqrt(2) - 1. We approximate sqrt(2) with 17/12
    tan_approx_num = sqrt2_num - sqrt2_den
    tan_approx_den = sqrt2_den
    print(f"The angle is 1/4 of a right angle. For this angle, tan(alpha) = sqrt(2) - 1.")
    integers_used.extend([1, 4, 2])
    print(f"To keep numbers simple, we can approximate sqrt(2) as the fraction {sqrt2_num}/{sqrt2_den}.")
    integers_used.extend([sqrt2_num, sqrt2_den])
    print(f"So, tan(alpha) is approx. {sqrt2_num}/{sqrt2_den} - 1 = {tan_approx_num}/{tan_approx_den}.")
    integers_used.extend([1, tan_approx_num, tan_approx_den])

    # h = d * tan(alpha)
    height = int(distance * tan_approx_num / tan_approx_den)
    print(f"The height h = {distance} * {tan_approx_num}/{tan_approx_den} = {int(distance/tan_approx_den)} * {tan_approx_num} = {height}m.")
    integers_used.extend([int(distance/tan_approx_den), height])
    print("-" * 30)

    # --- Part 2: Calculate the time (t) ---
    print("Part 2: Calculate the fall time (t)")
    print("The formula for fall time is t = sqrt(2 * h / g).")
    integers_used.append(2)
    print(f"We use our calculated height h = {height}m.")
    print(f"We use a simple value for gravity, g = {g} m/s^2.")
    integers_used.append(g)

    # t^2 = 2 * h / g
    t_squared = int(2 * height / g)
    print(f"So, t^2 = (2 * {height}) / {g} = {2 * height} / {g} = {t_squared}.")
    integers_used.append(2 * height)
    integers_used.append(t_squared)

    # t = sqrt(t_squared)
    print(f"This means t = sqrt({t_squared}). To simplify, we can write sqrt({t_squared}) = sqrt(4 * 5) = 2 * sqrt(5).")
    integers_used.extend([4, 5])
    
    # We approximate sqrt(5) with 9/4
    print(f"Now, we approximate sqrt(5) with the fraction {sqrt5_num}/{sqrt5_den}.")
    integers_used.extend([sqrt5_num, sqrt5_den])
    
    final_t_num = 2 * sqrt5_num
    final_t_den = sqrt5_den
    simplified_t_num = int(final_t_num / 2)
    simplified_t_den = int(final_t_den / 2)
    
    print(f"So, t is approx. 2 * {sqrt5_num}/{sqrt5_den} = {final_t_num}/{final_t_den}.")
    integers_used.extend([final_t_num, final_t_den])
    
    print(f"This fraction simplifies to {simplified_t_num}/{simplified_t_den}.")
    print(f"The final calculated time is {simplified_t_num / simplified_t_den} seconds.")
    print("-" * 30)

    # --- Final verification and answer ---
    # This part is for verification and not shown to the user
    # true_h = 240 * (math.sqrt(2) - 1)
    # true_t = math.sqrt(2 * true_h / 9.8)
    # approx_t = simplified_t_num / simplified_t_den
    # error = abs(true_t - approx_t)
    # print(f"DEBUG: Error is {error}, which should be < 0.1")

    largest_integer = max(integers_used)
    # The final answer format is Yz, where z is the largest integer.
    # This is not printed by the script as per instructions.
    # The final response will contain <<<Y{largest_integer}>>>

calculate_and_explain()