import numpy as np

def calculate_l(d, lambda_val):
    """
    Calculates the value of l(d, lambda) based on the derived formula.
    """

    # Define the function g(t) which is proportional to the unnormalized density of v
    # g(t) = exp(-t^2 / (2λ)) * (sin(t) / t)^(d-2)
    # We use np.sinc(t/np.pi) which is equivalent to sin(t)/t and handles t=0 correctly.
    def g(t, d, lambda_val):
        if t == 0:
            sinc_val = 1.0
        else:
            # np.sinc(x) is sin(pi*x)/(pi*x), so we use x = t/pi
            sinc_val = np.sinc(t / np.pi)
        
        term1 = np.exp(-t**2 / (2 * lambda_val))
        term2 = sinc_val**(d - 2)
        return term1 * term2

    # For x1, the norms of the pre-images v1 are t_1a and t_1b
    t_1a = np.pi / 2
    # Ensure the argument for arccos is valid
    arg1 = np.sqrt(3 / d)
    if arg1 > 1.0:
        print(f"Error: Invalid argument for arccos with d={d}. sqrt(3/d) > 1.")
        return None
    t_1b = np.arccos(arg1)

    # For x2, the norms of the pre-images v2 are t_2a and t_2b
    t_2a = np.pi / 2 # This is the same as t_1a
    # Ensure the argument for arccos is valid
    arg2 = np.sqrt(2 / d)
    if arg2 > 1.0:
        print(f"Error: Invalid argument for arccos with d={d}. sqrt(2/d) > 1.")
        return None
    t_2b = np.arccos(arg2)

    # Calculate the g values for each norm
    g_common = g(t_1a, d, lambda_val)
    g_1 = g(t_1b, d, lambda_val)
    g_2 = g(t_2b, d, lambda_val)

    # Calculate the numerator and denominator of the ratio
    numerator = g_common + g_1
    denominator = g_common + g_2
    
    # Calculate the final value l(d, lambda)
    if numerator <= 0 or denominator <= 0:
        print("Error: Numerator or denominator is non-positive, cannot take log.")
        return None
    ratio = numerator / denominator
    l_value = np.log(ratio)

    # Print the numbers in the final equation
    print(f"Calculating l(d, λ) for d = {d} and λ = {lambda_val}")
    print("-" * 30)
    print(f"The equation is: l = ln( (A + B) / (A + C) )")
    print(f"A = g(π/2) = {g_common}")
    print(f"B = g(arccos(sqrt(3/{d}))) = g({t_1b}) = {g_1}")
    print(f"C = g(arccos(sqrt(2/{d}))) = g({t_2b}) = {g_2}")
    print(f"Numerator (A + B) = {numerator}")
    print(f"Denominator (A + C) = {denominator}")
    print(f"Ratio = {ratio}")
    print(f"Final value l({d}, {lambda_val}) = {l_value}")
    
    return l_value

# Since d and lambda are not specified, we use example values d=4, lambda=1.
d_val = 4
lambda_val = 1
final_answer = calculate_l(d_val, lambda_val)

# The final answer is printed by the function above.
# For the required format, we will output the final numerical value here.
# print(f"\n<<< {final_answer} >>>")