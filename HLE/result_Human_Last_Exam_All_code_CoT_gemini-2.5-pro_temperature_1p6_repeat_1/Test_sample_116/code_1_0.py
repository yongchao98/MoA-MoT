import math

def solve_for_c(d, p):
    """
    Calculates the smallest constant c for the given inequality.

    The problem is to find the smallest constant c such that for any matrix
    A in R^(n x d) and any vector x in R^d, the following inequality holds:
    ||W^(1/2 - 1/p) * A*x||_2 <= c * ||A*x||_p
    where W is the diagonal matrix of the Lp Lewis weights of A.

    The solution from advanced matrix theory is c = d**|1/2 - 1/p|.
    This function computes and prints this value.

    Args:
        d (int): The number of columns of the matrix A (a fixed integer > 0).
        p (float): The parameter for the Lp norm (a float in (0, inf)).
    """
    if d <= 0 or not isinstance(d, int):
        print("Error: d must be a positive integer.")
        return
    if p <= 0:
        print("Error: p must be a positive number.")
        return

    exponent = abs(0.5 - 1/p)
    c = d ** exponent

    # Print the equation with the numbers substituted in.
    print(f"For d = {d} and p = {p}:")
    print(f"The constant c is calculated by the formula: d^(|1/2 - 1/p|)")
    print(f"c = {d}^(|0.5 - 1/{p}|)")
    print(f"c = {d}^(|0.5 - {1/p:.4f}|)")
    print(f"c = {d}^({exponent:.4f})")
    print(f"c = {c:.4f}")

# Example usage with d=10 and p=4, as discussed in the thought process.
# You can change these values to test other cases.
d_val = 10
p_val = 4.0
solve_for_c(d_val, p_val)

print("\n" + "-"*20 + "\n")

# Another example with p < 2.
d_val = 10
p_val = 1.5
solve_for_c(d_val, p_val)
