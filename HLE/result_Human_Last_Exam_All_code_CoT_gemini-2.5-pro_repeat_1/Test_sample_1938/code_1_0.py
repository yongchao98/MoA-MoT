import math

def solve():
    """
    Calculates the minimum q_0 such that the probability of Alice winning is at least 95%,
    and returns floor(100 * q_0).
    """
    # The target winning probability for Alice.
    p = 0.95

    # The equation to solve for q_0 is:
    # q_0^4 = (1 - (1 - p)^(1/3)) / p^3
    
    one_minus_p = 1 - p

    # Calculate the numerator and denominator of the fraction.
    numerator = 1 - (one_minus_p)**(1/3)
    denominator = p**3

    # Calculate q_0^4 and then q_0.
    q_pow_4 = numerator / denominator
    q0 = q_pow_4**(1/4)

    # The final result is floor(100 * q_0).
    result = math.floor(100 * q0)
    
    print("Step-by-step calculation:")
    print(f"The equation for q_0 is: q_0 = ((1 - (1 - p)^(1/3)) / p^3)^(1/4)")
    print(f"Given p = {p}, we have:")
    print(f"  1 - p = {one_minus_p}")
    print(f"  (1 - p)^(1/3) = {one_minus_p**(1/3):.6f}")
    print(f"  Numerator = 1 - {one_minus_p**(1/3):.6f} = {numerator:.6f}")
    print(f"  p^3 = {denominator:.6f}")
    print(f"  q_0^4 = {numerator:.6f} / {denominator:.6f} = {q_pow_4:.6f}")
    print(f"  q_0 = ({q_pow_4:.6f})^(1/4) = {q0:.6f}")
    print("-" * 20)
    print(f"The value of 100 * q_0 is: 100 * {q0:.6f} = {100*q0:.6f}")
    print(f"The final answer is floor(100 * q_0) = {result}")

solve()