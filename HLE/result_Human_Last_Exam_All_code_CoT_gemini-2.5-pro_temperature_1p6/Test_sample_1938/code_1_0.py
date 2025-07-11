import math

def solve_for_q():
    """
    Calculates the minimum q for Alice to win with at least 95% probability.
    """
    # The problem boils down to solving the equation x = 1 - (1 - q^4 * x^3)^3
    # for q, given that Alice's winning probability x is 0.95.

    # Alice's winning probability
    x = 0.95
    print(f"We are given Alice's winning probability, x = {x}\n")
    
    # We rearrange the equation to solve for q:
    # 1 - x = (1 - q^4 * x^3)^3
    # (1 - x)^(1/3) = 1 - q^4 * x^3
    # q^4 * x^3 = 1 - (1 - x)^(1/3)
    # q^4 = (1 - (1 - x)^(1/3)) / x^3
    # q = [ (1 - (1 - x)^(1/3)) / x^3 ]^(1/4)
    print("The equation to find q is: q = [(1 - (1 - x)^(1/3)) / x^3]^(1/4)\n")

    print("--- Calculating components of the equation ---")
    
    # Calculate 1 - x
    one_minus_x = 1 - x
    print(f"Value of (1 - x): {one_minus_x:.4f}")
    
    # Calculate (1 - x)^(1/3)
    cbrt_one_minus_x = one_minus_x**(1/3)
    print(f"Value of (1 - x)^(1/3): {cbrt_one_minus_x:.4f}")
    
    # This is the numerator of our main fraction
    numerator = 1 - cbrt_one_minus_x
    print(f"Value of the numerator [1 - (1 - x)^(1/3)]: {numerator:.4f}")

    # This is the denominator of our main fraction
    x_cubed = x**3
    print(f"Value of the denominator [x^3]: {x_cubed:.4f}")
    
    # Calculate q^4
    q_pow_4 = numerator / x_cubed
    print(f"\nValue of q^4: {q_pow_4:.4f}")

    # Finally, calculate q
    q_0 = q_pow_4**(1/4)
    print(f"The minimum required value of q (q_0) is: {q_0:.4f}")
    
    # The question asks for floor(100 * q_0)
    final_answer = math.floor(100 * q_0)
    print(f"\nThe value of 100 * q_0 is: {100 * q_0:.4f}")
    print(f"The final answer, floor(100 * q_0), is: {final_answer}")

solve_for_q()