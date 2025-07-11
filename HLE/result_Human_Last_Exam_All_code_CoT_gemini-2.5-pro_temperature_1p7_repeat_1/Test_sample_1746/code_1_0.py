def check_division_identity(n, m):
    """
    Checks the Division Algorithm Identity for integers n and m.
    The identity is: n == (n // m) * m + (n % m)
    """
    print(f"Checking for n = {n}, m = {m}")

    quotient = n // m
    remainder = n % m
    
    # Calculate the right-hand side of the equation
    rhs = quotient * m + remainder
    
    # Output each number in the final equation
    print(f"({n} // {m}) * {m} + ({n} % {m}) == {n}")
    print(f"({quotient}) * {m} + ({remainder}) == {n}")
    print(f"{quotient * m} + {remainder} == {n}")
    print(f"{rhs} == {n}")
    
    # Check if the identity holds true
    is_true = (rhs == n)
    print(f"Result: {is_true}\n")

# Per statement G, check for positive integers
check_division_identity(17, 5)

# The identity also holds true for negative integers, showcasing its robustness
check_division_identity(-17, 5)
check_division_identity(17, -5)
check_division_identity(-17, -5)
