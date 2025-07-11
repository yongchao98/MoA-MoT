import math

def solve():
    """
    This function calculates the minimum q such that Alice's winning probability is at least 95%.
    """
    
    # Let p_A be the probability of Alice winning. We are given p_A >= 0.95.
    # The minimum value of q (q_0) will occur when p_A is exactly 0.95.
    p_A = 0.95
    
    # The relationship between p_A and q is derived from the game's mechanics:
    # p_A = 1 - (1 - q^4 * p_A^3)^3
    # We need to rearrange this equation to solve for q.
    
    # 1 - p_A = (1 - q^4 * p_A^3)^3
    # (1 - p_A)^(1/3) = 1 - q^4 * p_A^3
    # q^4 * p_A^3 = 1 - (1 - p_A)^(1/3)
    # q^4 = (1 - (1 - p_A)^(1/3)) / p_A^3
    # q = ((1 - (1 - p_A)^(1/3)) / p_A^3)^(1/4)
    
    print("Let p_A be the probability Alice wins. We are solving for q when p_A = 0.95.")
    print("The equation to solve is: p_A = 1 - (1 - q^4 * p_A^3)^3")
    print(f"Rearranging to solve for q, we get: q = ((1 - (1 - p_A)^(1/3)) / p_A^3)^(1/4)")
    print("-" * 30)
    print("Calculation steps:")

    pa_cubed = p_A**3
    print(f"p_A^3 = ({p_A})^3 = {pa_cubed}")
    
    one_minus_pa = 1 - p_A
    print(f"1 - p_A = 1 - {p_A} = {one_minus_pa}")
    
    cbrt_one_minus_pa = one_minus_pa**(1/3)
    print(f"(1 - p_A)^(1/3) = ({one_minus_pa})^(1/3) = {cbrt_one_minus_pa}")
    
    numerator = 1 - cbrt_one_minus_pa
    print(f"Numerator: 1 - {cbrt_one_minus_pa} = {numerator}")

    q_pow_4 = numerator / pa_cubed
    print(f"q^4 = Numerator / p_A^3 = {numerator} / {pa_cubed} = {q_pow_4}")
    
    q0 = q_pow_4**(1/4)
    print(f"q = (q^4)^(1/4) = ({q_pow_4})^(1/4) = {q0}")
    print("-" * 30)

    # The result required is floor(100 * q_0)
    final_value = 100 * q0
    result = math.floor(final_value)
    
    print(f"The minimum value of q is q_0 ≈ {q0:.5f}")
    print(f"We need to calculate floor(100 * q_0).")
    print(f"100 * q_0 ≈ {final_value:.5f}")
    print(f"The final answer is floor({final_value:.5f}) = {result}")

solve()