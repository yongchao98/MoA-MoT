import math

def explain_calculation():
    """
    This function prints a step-by-step explanation of the calculation
    in a way a child can understand, using small integers and fractions.
    """
    
    # Initial known values
    distance = 240
    
    # Approximations for simple calculation
    g_approx = 10
    tan_alpha_num = 5
    tan_alpha_den = 12

    # Break down distance to keep numbers small
    d_mantissa = 24
    d_exponent = 10

    # Start of the explanation
    print("Yes, your son can solve this! Here is a simple way to do the calculation.")
    print("-" * 30)
    print("First, we need the formula to find the time squared (t * t).")
    print(f"The formula is: t * t = (2 * distance * tan(alpha)) / g")
    print(f"\nWe will use simple numbers for the calculation:")
    print(f"The distance is {distance}m. Let's write it as {d_mantissa} * {d_exponent} to make it easier.")
    print(f"For gravity, g, we can use the number {g_approx}.")
    print(f"For the angle part, tan(alpha), we can use the fraction {tan_alpha_num}/{tan_alpha_den}.")
    
    print("\nNow let's put the numbers in the formula:")
    
    # This is the key step, showing all the numbers in the equation
    print(f"t * t = (2 * ({d_mantissa} * {d_exponent}) * {tan_alpha_num}/{tan_alpha_den}) / {g_approx}")
    
    print("\nLook! We are multiplying by 10 and dividing by 10, so they cancel each other out.")
    print(f"The calculation becomes much simpler:")
    
    final_eq_part1 = 2
    final_eq_part2 = d_mantissa
    final_eq_part3 = tan_alpha_num
    final_eq_part4 = tan_alpha_den
    
    print(f"t * t = ({final_eq_part1} * {final_eq_part2} * {final_eq_part3}) / {final_eq_part4}")
    
    intermediate_result1 = final_eq_part2 // final_eq_part4
    print(f"\nFirst, we can calculate {final_eq_part2} / {final_eq_part4}, which is {intermediate_result1}.")
    
    final_result_sq = final_eq_part1 * intermediate_result1 * final_eq_part3
    print(f"So now we have: t * t = {final_eq_part1} * {intermediate_result1} * {final_eq_part3}")
    print(f"This means: t * t = {final_result_sq}")
    
    print("\nNow we need to find the time, t.")
    print(f"We know 4 * 4 = 16 and 5 * 5 = 25. The number we are looking for is between 4 and 5.")
    print("A good guess is 4 and a half, which is the fraction 9/2.")
    
    t_approx_num = 9
    t_approx_den = 2
    check_num = t_approx_num * t_approx_num
    check_den = t_approx_den * t_approx_den
    
    print(f"Let's check: ({t_approx_num}/{t_approx_den}) * ({t_approx_num}/{t_approx_den}) = {check_num}/{check_den}")
    print(f"{check_num}/{check_den} is 20.25, which is very close to {final_result_sq}!")
    
    print(f"\nSo, the time for the egg to drop is about {t_approx_num}/{t_approx_den} seconds.")

explain_calculation()