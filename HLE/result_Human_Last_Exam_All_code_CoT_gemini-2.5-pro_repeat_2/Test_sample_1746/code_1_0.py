def demonstrate_division_identity():
    """
    This function demonstrates the truth of statement G:
    "For positive integers n, m: (n//m)*m + n%m == n holds true in all versions"

    This identity is a fundamental property of Euclidean division for integers.
    Python's floor division (//) and modulo (%) operators are defined to
    satisfy this relationship for all integers (and in this case, we test positive ones).
    This has been consistent across all Python versions.
    """
    n = 27
    m = 4

    # Calculate the quotient and remainder
    quotient = n // m
    remainder = n % m

    # Calculate the left-hand side of the equation
    lhs_result = quotient * m + remainder

    print("Demonstrating the identity: (n // m) * m + n % m == n")
    print(f"For n = {n} and m = {m}:\n")

    # Print the equation with variables
    print(f"({n} // {m}) * {m} + {n} % {m} == {n}")

    # Print the equation with calculated intermediate values
    print(f"({quotient}) * {m} + {remainder} == {n}")

    # Print the final calculated result and the boolean outcome
    print(f"{lhs_result} == {n}, which is {lhs_result == n}")

demonstrate_division_identity()