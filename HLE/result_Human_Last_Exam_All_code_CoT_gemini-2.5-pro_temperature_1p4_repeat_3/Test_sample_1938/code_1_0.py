import math

def solve():
    """
    Calculates the minimum q for Alice's win probability to be at least 95%.
    """
    # The target probability for Alice to win.
    p_A = 0.95
    
    # The equation relating p_A and q is p_A = 1 - (1 - q^4 * p_A^3)^3.
    # We solve for q_0 when p_A = 0.95.
    # q_0 = ((1 - (1 - p_A)^(1/3)) / p_A^3)^(1/4)
    
    p_A_cubed = p_A**3
    one_minus_p_A = 1 - p_A
    cbrt_one_minus_p_A = one_minus_p_A**(1/3)
    
    numerator = 1 - cbrt_one_minus_p_A
    denominator = p_A_cubed
    
    q0_pow_4 = numerator / denominator
    q0 = q0_pow_4**(1/4)
    
    # The final answer is floor(100 * q0)
    result = math.floor(100 * q0)
    
    print("The final equation to solve for q is: q = ((1 - (1 - p)^(1/3)) / p^3)^(1/4)")
    print("We substitute p = 0.95 into this equation.")
    print(f"p = {p_A}")
    print(f"p^3 = {denominator}")
    print(f"1 - p = {one_minus_p_A}")
    print(f"(1 - p)^(1/3) = {cbrt_one_minus_p_A}")
    print(f"1 - (1 - p)^(1/3) = {numerator}")
    print(f"q^4 = {numerator} / {denominator} = {q0_pow_4}")
    print(f"q = ({q0_pow_4})^(1/4) = {q0}")
    print(f"The required value is floor(100 * q).")
    print(f"100 * q = {100 * q0}")
    print(f"floor(100 * q) = {result}")

solve()
