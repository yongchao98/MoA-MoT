import math

def calculate_l(n: int):
    """
    Calculates the exact value of the function l(n) for n >= 5.

    The function is defined as l(n) = f^(1)(f_M^(2)(f^(3)(P))), where the
    intermediate calculations lead to the final formula implemented below.
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("Input n must be an integer greater than or equal to 5.")

    # The calculation of l(n) simplifies to:
    # l(n) = Sum of first and last rows of Q - Sum of first and last rows of N
    # where Q is the inverse of P, and N is the projection of Q onto the normal space.
    
    # The sum of the first and last rows of Q is 6.
    sum_rows_Q = 6.0

    # The parameters 'a' and 'b' from the matrix M_n
    a = math.sqrt(1 - (n - 1) / n**2)
    b = 1 / n

    # The sum of the first and last rows of the normal projection component N
    # simplifies to the following expression:
    # 4 * (a^2 + (n - 1/2)ab + (n - 3/2)b^2)
    term_a_sq = a**2
    term_ab = (n - 0.5) * a * b
    term_b_sq = (n - 1.5) * b**2
    
    sum_rows_N = 4 * (term_a_sq + term_ab + term_b_sq)

    # The final value of l(n)
    l_n = sum_rows_Q - sum_rows_N
    
    # The instruction "output each number in the final equation" is interpreted
    # as showing the main components of the calculation.
    print(f"For n = {n}:")
    print(f"Sum of relevant rows of Q = {sum_rows_Q}")
    print(f"Sum of relevant rows of the projection term N = {sum_rows_N}")
    print(f"l({n}) = {sum_rows_Q} - {sum_rows_N} = {l_n}")

# Example calculation for n=5
calculate_l(5)

# Example calculation for n=10
calculate_l(10)
