import numpy as np

def calculate_bound(n):
    """
    Calculates the upper bound of the product for a given n.
    The bound is given by the formula 2 * (n + 2) / (n + 1).
    """
    if n <= 0:
        raise ValueError("n must be a positive integer.")
    
    numerator = n + 2
    denominator = n + 1
    bound = 2 * numerator / denominator
    return bound

def main():
    """
    Calculates the least upper bound of the product over all positive integers n.
    """
    print("Under the assumption that the problem refers to an orthogonal tridiagonalization,")
    print("the product is bounded by the function f(n) = 2 * (n+2) / (n+1).")
    print("\nCalculating the bound for the first 10 positive integers:")

    for n_val in range(1, 11):
        bound_val = calculate_bound(n_val)
        print(f"For n = {n_val:2d}, the upper bound is 2 * ({n_val + 2}/{n_val + 1}) = {bound_val:.4f}")

    # The function f(n) = 2 * (1 + 1/(n+1)) is a decreasing function for n > 0.
    # The least upper bound (supremum) is its value at n=1.
    n_lub = 1
    lub_value = calculate_bound(n_lub)

    print("\nAs the sequence of bounds is decreasing, the least upper bound is the value at n=1.")
    print(f"LUB = 2 * ({n_lub} + 2) / ({n_lub} + 1) = {int(lub_value)}")

if __name__ == "__main__":
    main()
