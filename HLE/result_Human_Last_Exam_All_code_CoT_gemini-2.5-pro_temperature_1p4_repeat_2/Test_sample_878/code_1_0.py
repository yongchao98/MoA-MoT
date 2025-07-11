import math

def calculate_egg_drop_time():
    """
    This function calculates the time for an egg to fall from a skyscraper,
    using approximations suitable for a young student.
    It prints the calculation step-by-step.
    """

    # Given values from the problem
    distance = 240  # meters

    # Approximations for the calculation
    # g is approximated as 10 m/s^2 for simplicity.
    g = 10
    # tan(22.5 degrees) is approximated as the fraction 27/64.
    # The true value is sqrt(2)-1 which is approx 0.414. 27/64 is approx 0.422.
    tan_alpha_num = 27
    tan_alpha_den = 64

    print("Yes, your son can solve this. Here is a possible calculation with small integers:")
    print("-" * 30)
    print("The formula for the time 't' to fall from a height 'h' is: t^2 = (2 * h) / g")
    print("The height 'h' is calculated as: h = distance * tan(alpha)")
    print(f"The given distance is {distance}m and we will use approximations for g = {g} and tan(alpha) = {tan_alpha_num}/{tan_alpha_den}.")
    print("-" * 30)

    # Let's write the full expression for t^2
    print("1. First, let's write the full equation for t^2:")
    print(f"t^2 = (2 * distance * tan(alpha)) / g")
    print(f"t^2 = (2 * {distance} * ({tan_alpha_num} / {tan_alpha_den})) / {g}")
    print(f"t^2 = (2 * {distance} * {tan_alpha_num}) / ({g} * {tan_alpha_den})\n")

    # To keep the numbers small, we can simplify the expression step by step
    # instead of multiplying everything at once.
    print("2. To make it easier, let's simplify part of the fraction first:")
    term1_num = 2 * distance
    term1_den = g
    term1_val = term1_num // term1_den
    print(f"   (2 * {distance}) / {g}  = {term1_num} / {term1_den} = {term1_val}")
    print(f"So, the equation becomes: t^2 = {term1_val} * ({tan_alpha_num} / {tan_alpha_den})\n")

    print("3. Now, let's write this as a single fraction:")
    print(f"t^2 = ({term1_val} * {tan_alpha_num}) / {tan_alpha_den}")
    # Instead of calculating 48 * 27, we can simplify the fraction 48/64 first.
    common_divisor = math.gcd(term1_val, tan_alpha_den)
    simplified_num1 = term1_val // common_divisor
    simplified_den1 = tan_alpha_den // common_divisor
    print(f"   We can simplify the fraction by dividing {term1_val} and {tan_alpha_den} by their greatest common divisor, which is {common_divisor}.")
    print(f"   {term1_val} / {common_divisor} = {simplified_num1}")
    print(f"   {tan_alpha_den} / {common_divisor} = {simplified_den1}")
    print(f"So, the equation is now: t^2 = ({simplified_num1} / {simplified_den1}) * {tan_alpha_num}\n")

    print("4. Now we can multiply the remaining numbers:")
    final_num = simplified_num1 * tan_alpha_num
    final_den = simplified_den1
    print(f"t^2 = ({simplified_num1} * {tan_alpha_num}) / {simplified_den1}")
    print(f"t^2 = {final_num} / {final_den}\n")
    
    print("5. The final step is to take the square root:")
    # The numerator and denominator are perfect squares
    t_num = int(math.sqrt(final_num))
    t_den = int(math.sqrt(final_den))
    print(f"t = sqrt({final_num}) / sqrt({final_den})")
    print(f"t = {t_num} / {t_den}\n")

    final_time = t_num / t_den
    print("-" * 30)
    print(f"The calculated time is {final_time} seconds.")
    print("-" * 30)

calculate_egg_drop_time()