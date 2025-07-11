import sys

def verify_division_identity():
    """
    This function verifies statement G, which claims that for positive integers
    n and m, the identity (n//m)*m + n%m == n holds true. This is a
    fundamental property of integer division and modulo in Python.
    """
    n = 27
    m = 5
    
    # Perform the floor division and modulo operations
    quotient = n // m
    remainder = n % m
    
    # Calculate the left-hand side of the equation
    lhs_result = quotient * m + remainder
    
    # Print the equation step by step as requested
    print("Verifying the identity: (n//m)*m + n%m == n")
    print(f"For n = {n} and m = {m}:\n")
    
    # Print the first line of the equation with the operators
    print(f"({n}//{m})*{m} + {n}%{m} == {n}")
    
    # Print the second line with the results of the division and modulo
    print(f"({quotient})*{m} + {remainder} == {n}")

    # Print the final numerical comparison
    print(f"{lhs_result} == {n}")
    
    # Final confirmation
    is_true = lhs_result == n
    print(f"\nThe statement (n//m)*m + n%m == n is: {is_true}")

verify_division_identity()