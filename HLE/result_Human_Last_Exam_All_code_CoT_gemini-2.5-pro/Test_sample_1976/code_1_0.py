import math

def calculate_T_norm_1(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n with odd n.

    Args:
        n (int): An odd integer.

    Returns:
        float: The 1-norm of the correlation matrix T.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")

    n_plus_1 = n + 1
    denominator = 2**n_plus_1
    result = 1 - 1 / denominator

    print(f"For n = {n}:")
    print(f"The formula for the 1-norm of T is 1 - 1 / 2^(n+1).")
    print(f"Substituting n = {n}, we get:")
    print(f"1 - 1 / 2^({n}+1) = 1 - 1 / 2^{n_plus_1} = 1 - 1 / {denominator} = {(denominator-1)}/{denominator} = {result}")
    print("-" * 20)
    return result

if __name__ == '__main__':
    # Demonstrate the calculation for a few odd values of n.
    odd_n_values = [1, 3, 5, 13]
    for n_val in odd_n_values:
        calculate_T_norm_1(n_val)
