import math

def compute_c(d, p):
    """
    Computes the smallest constant c for the given inequality.

    Args:
        d (int): The dimension of the subspace, d > 0.
        p (float): The L_p norm exponent, p > 0.

    Returns:
        float: The value of the constant c.
    """
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(p, (int, float)) or p <= 0:
        raise ValueError("p must be a positive number.")

    if p >= 2:
        # For p >= 2, c = d^(1/2 - 1/p)
        val = d**(0.5 - 1/p)
    else:
        # For 0 < p < 2, c = 1
        val = 1.0
    
    return val

def main():
    """
    Main function to demonstrate the calculation of c.
    """
    # Example 1: d=4, p=4
    d1 = 4
    p1 = 4
    c1 = compute_c(d1, p1)
    print(f"For d = {d1} and p = {p1}:")
    print(f"The formula is c = d**(1/2 - 1/p) = {d1}**(1/2 - 1/{p1})")
    print(f"The constant c is: {c1}")
    print("-" * 20)

    # Example 2: d=10, p=1
    d2 = 10
    p2 = 1
    c2 = compute_c(d2, p2)
    print(f"For d = {d2} and p = {p2}:")
    print(f"The formula is c = 1")
    print(f"The constant c is: {c2}")
    print("-" * 20)
    
    # Example 3: d=9, p=2
    d3 = 9
    p3 = 2
    c3 = compute_c(d3, p3)
    print(f"For d = {d3} and p = {p3}:")
    # The two formulas coincide here
    print(f"The formula gives c = max(1, {d3}**(1/2 - 1/{p3})) = max(1, 1)")
    print(f"The constant c is: {c3}")
    print("-" * 20)

if __name__ == '__main__':
    main()