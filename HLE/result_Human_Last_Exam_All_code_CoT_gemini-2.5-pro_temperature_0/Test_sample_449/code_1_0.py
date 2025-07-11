import math
import numpy

def solve_probability():
    """
    Calculates the probability that the conditioned random walk avoids the neighbors of the origin.
    """
    # --- Step 1: Define constants and initial parameters ---
    pi = math.pi
    # Euler-Mascheroni constant
    gamma = numpy.euler_gamma
    start_point = (3000, 4000)
    
    # --- Step 2: Calculate the norm of the starting vector ---
    norm_x0 = math.sqrt(start_point[0]**2 + start_point[1]**2)
    
    # --- Step 3: Calculate the potential kernel h(x) for the set N ---
    # The set N are the neighbors of the origin, all at norm 1.
    # The potential kernel a(1,0) has a known exact value.
    h_N = 4 / pi - 1
    
    # --- Step 4: Calculate the potential kernel h(x) for the starting point x0 ---
    # For large ||x||, we use the asymptotic formula:
    # a(x) ~ (2/pi) * (ln(||x||) + gamma + ln(8))
    log_term = math.log(norm_x0)
    constant_term = gamma + math.log(8)
    h_x0 = (2 / pi) * (log_term + constant_term)
    
    # --- Step 5: Calculate the final escape probability ---
    # The probability is given by the formula P = 1 - h_N / h_x0
    probability = 1 - h_N / h_x0
    
    # --- Step 6: Print the results as requested ---
    print("The problem is to calculate the escape probability P = 1 - h_N / h_x0, where:")
    print(f" - h_N is the potential kernel value on the set of the origin's neighbors.")
    print(f" - h_x0 is the potential kernel value at the starting point ({start_point[0]},{start_point[1]}).")
    print("\nCalculating the values:")
    print(f"h_N = a(1,0) = 4/pi - 1 = {h_N:.5f}")
    print(f"h_x0 = a({start_point[0]},{start_point[1]}) ≈ (2/pi) * (ln({norm_x0}) + gamma + ln(8)) = {h_x0:.5f}")
    
    print("\nThe final equation with these numbers is:")
    print(f"P = 1 - {h_N:.5f} / {h_x0:.5f}")
    
    print("\nFinal result:")
    print(f"P = {probability:.5f}")
    
    # Approximate answer with two significant digits
    approx_answer = float(f"{probability:.2g}")
    print(f"\nThe approximate answer with two significant digits is: {approx_answer}")
    
    return approx_answer

if __name__ == '__main__':
    solve_probability()
