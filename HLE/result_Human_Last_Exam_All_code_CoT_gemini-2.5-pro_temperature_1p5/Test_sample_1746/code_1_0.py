def verify_division_identity():
    """
    This function verifies statement G:
    For positive integers n, m: (n//m)*m + n%m == n holds true in all versions.
    
    This identity is the mathematical definition of how Python's floor division
    and modulo operators relate to each other. It is guaranteed by the language
    to hold true for all integers.
    """
    print("Verifying statement G: (n//m)*m + n%m == n")
    print("-" * 50)
    
    # A set of test cases with positive integers
    test_cases = [(10, 3), (20, 7), (42, 6), (5, 8), (12345, 67)]
    
    all_passed = True
    
    for n, m in test_cases:
        print(f"Testing with n = {n}, m = {m}")
        
        # Calculate each component of the equation
        quotient = n // m
        remainder = n % m
        left_hand_side = quotient * m + remainder
        right_hand_side = n
        
        # Check if the identity holds
        result = (left_hand_side == right_hand_side)
        
        if not result:
            all_passed = False
            
        # As requested, output each number in the final equation step-by-step
        print(f"Equation: ({n} // {m}) * {m} + ({n} % {m}) == {n}")
        print(f"Step 1:   ({quotient}) * {m} + ({remainder}) == {n}")
        print(f"Step 2:   {quotient * m} + {remainder} == {n}")
        print(f"Result:   {left_hand_side} == {right_hand_side} -> {result}")
        print("-" * 20)
        
    print(f"\nConclusion: The statement held true for all test cases: {all_passed}")

verify_division_identity()