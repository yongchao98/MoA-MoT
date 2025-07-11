def demonstrate_division_identity():
    """
    This function demonstrates the truth of statement (G):
    For positive integers n, m: (n//m)*m + n%m == n holds true in all versions.

    This is the Division Algorithm identity, which defines the relationship
    between floor division ('//') and modulo ('%') operators for integers.
    """
    print("Verifying statement G: (n//m) * m + n % m == n")
    print("-" * 50)

    test_cases = [(17, 5), (42, 7), (100, 3)]

    for n, m in test_cases:
        quotient = n // m
        remainder = n % m
        result = quotient * m + remainder

        print(f"Testing with n = {n}, m = {m}:")

        # Output each number in the final equation step-by-step
        # Final equation format: quotient * m + remainder = result
        print(f"  ({n} // {m}) * {m} + ({n} % {m}) = {n}")
        print(f"  ({quotient}) * {m} + ({remainder}) = {n}")
        print(f"  {quotient * m} + {remainder} = {n}")
        print(f"  {result} = {n}")

        # Final check
        if result == n:
            print("  Conclusion: The identity holds true.\n")
        else:
            print("  Conclusion: The identity FAILED.\n")

demonstrate_division_identity()