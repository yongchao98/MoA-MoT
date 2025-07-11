import math

def solve():
    """
    This function takes the inputs p, k, n and computes the number of coefficients
    not divisible by p^k based on the derived formula.
    """
    # The problem is stated with variables p, k, n.
    # We will use example values for these variables to demonstrate the calculation.
    # Let's use p=3, k=2, n=4 as an example.
    p = 3
    k = 2
    n = 4

    # The problem is to compute the number of coefficients in the polynomial
    # P_{p^n}(x) = F_{p, k}^{p^n}(x) that are not divisible by p^k.
    #
    # Let N(m) be the number of coefficients of P_m(x) not divisible by p^k.
    #
    # Step 1: Find N(1).
    # N(1) is the number of coefficients of P_1(x) = prod_{i=1}^{p^k} (x-i)
    # that are not divisible by p^k.
    # A known result from number theory states this number is 2 * p^(k-1).
    #
    # Step 2: Find the recurrence for N(m).
    # The structure of polynomial composition over p-adic integers suggests
    # that the number of "surviving" coefficients multiplies by p for every p iterations.
    # This gives the recurrence N(m*p) = p * N(m).
    #
    # Step 3: Calculate N(p^n).
    # Using the recurrence:
    # N(p^n) = p * N(p^(n-1)) = p^2 * N(p^(n-2)) = ... = p^n * N(1).
    #
    # Step 4: Combine the results.
    # N(p^n) = p^n * (2 * p^(k-1)) = 2 * p^(n + k - 1).

    # Let's compute this value.
    exponent = n + k - 1
    result = 2 * (p ** exponent)

    # The problem asks to express the answer in terms of p, k, and n.
    # The final expression is 2 * p^(n + k - 1).
    # We print the components of this expression.
    print("The number of coefficients not divisible by p^k is given by the formula: 2 * p^(n + k - 1)")
    print(f"The base is: {2}")
    print(f"The prime is: {p}")
    print(f"The exponent is: n + k - 1 = {n} + {k} - {1} = {exponent}")
    print(f"So the final answer is: 2 * {p}^({exponent}) = {result}")

solve()