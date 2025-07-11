import math

def calculate_T_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for odd n.

    Args:
        n: An odd integer.
    """
    if n % 2 == 0:
        print(f"The provided formula is for odd n. Received n={n}.")
        return

    # The formula for the 1-norm is 4**(n+1) - 2**(n+1)
    term1 = 4**(n + 1)
    term2 = 2**(n + 1)
    result = term1 - term2

    print(f"For n = {n}:")
    print(f"The calculation is 4^({n+1}) - 2^({n+1}).")
    # Using the requirement to output each number in the final equation
    print(f"{term1} - {term2} = {result}")
    print(f"The 1-norm of the correlation matrix T is {result}.\n")

if __name__ == '__main__':
    # Calculate the norm for a few small odd values of n.
    print("Calculating the 1-norm of the correlation matrix T for J_n for some odd n.")
    calculate_T_norm(1)
    calculate_T_norm(3)
    calculate_T_norm(5)