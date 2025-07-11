import math

def solve_for_q():
    """
    This function solves for the minimum value of q (q_0) such that the
    probability of Alice winning is at least 95%, and then calculates
    floor(100 * q_0).
    """

    # The probability of Alice winning, p_A, is given as 0.95.
    p_A = 0.95

    # We derived the following relationship between p_A and q:
    # p_A = 1 - (1 - q^4 * p_A^3)^3
    #
    # Rearranging to solve for q, we get:
    # q = ((1 - (1 - p_A)^(1/3)) / p_A^3)^(1/4)
    #
    # Now, we will plug in the value of p_A and calculate q.

    print("Step-by-step calculation for q_0:")
    print("---------------------------------")
    print(f"Given probability for Alice to win, p_A = {p_A}")

    # Calculate 1 - p_A
    one_minus_p_A = 1 - p_A
    print(f"1 - p_A = {one_minus_p_A:.2f}")

    # Calculate the cube root of (1 - p_A)
    cbrt_val = one_minus_p_A**(1/3)
    print(f"(1 - p_A)^(1/3) = {cbrt_val}")

    # Calculate the numerator of the fraction for q^4
    numerator = 1 - cbrt_val
    print(f"Numerator: 1 - (1 - p_A)^(1/3) = {numerator}")

    # Calculate the denominator of the fraction for q^4
    denominator = p_A**3
    print(f"Denominator: p_A^3 = {denominator}")

    # Calculate q^4
    q_pow_4 = numerator / denominator
    print(f"q_0^4 = Numerator / Denominator = {q_pow_4}")

    # Calculate q_0
    q_0 = q_pow_4**(1/4)
    print(f"q_0 = (q_0^4)^(1/4) = {q_0}")

    # Calculate 100 * q_0
    result_float = 100 * q_0
    print(f"\n100 * q_0 = {result_float}")

    # Calculate the floor of the result
    final_answer = math.floor(result_float)
    print(f"The final answer, floor(100 * q_0), is: {final_answer}")

# Run the solver function
solve_for_q()