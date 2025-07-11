def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    n = 4
    q = 2

    # The order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.
    terms = []
    for i in range(n):
        term = q**n - q**i
        terms.append(term)

    # Calculate the final order by multiplying the terms.
    order = 1
    for term in terms:
        order *= term

    # Format the output string for the equation.
    equation_str = " * ".join(map(str, terms))

    print("The highest possible order for E is the order of Out(D), which is isomorphic to GL(4, F_2).")
    print("The order of GL(4, F_2) is calculated as:")
    print(f"(2^4 - 2^0) * (2^4 - 2^1) * (2^4 - 2^2) * (2^4 - 2^3)")
    print("This evaluates to the following equation:")
    print(f"{equation_str} = {order}")

solve()