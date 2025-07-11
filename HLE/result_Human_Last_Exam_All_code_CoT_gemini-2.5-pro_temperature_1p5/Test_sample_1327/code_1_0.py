def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a two-part pattern.
    The main pattern is a_n = 2^(n+1) + b_k
    The secondary pattern for b_k is b_k = 3 * b_{k-1} - 2
    """
    sequence = [2, 11, 23, 51, 119]
    print(f"Given sequence: {sequence}\n")

    print("Step 1: Discovering the pattern.")
    print("The pattern (starting from the 2nd term) is: a_n = 2^(n+1) + b, where b follows its own pattern.")
    
    # Verify the pattern for existing terms
    # n=2: a_2 = 11. b is 11 - 2^(2+1) = 11 - 8 = 3
    b = 3
    print(f"For n=2: {sequence[1]} = 2^(2+1) + {b} = 8 + {b}")

    # n=3: a_3 = 23. Next b is 3*3 - 2 = 7.
    # Check: 2^(3+1) + 7 = 16 + 7 = 23. Correct.
    b_prev = b
    b = 3 * b_prev - 2
    print(f"For n=3: {sequence[2]} = 2^(3+1) + (3*{b_prev}-2) = 16 + {b}")

    # n=4: a_4 = 51. Next b is 3*7 - 2 = 19.
    # Check: 2^(4+1) + 19 = 32 + 19 = 51. Correct.
    b_prev = b
    b = 3 * b_prev - 2
    print(f"For n=4: {sequence[3]} = 2^(4+1) + (3*{b_prev}-2) = 32 + {b}")

    # n=5: a_5 = 119. Next b is 3*19 - 2 = 55.
    # Check: 2^(5+1) + 55 = 64 + 55 = 119. Correct.
    b_prev = b
    b = 3 * b_prev - 2
    print(f"For n=5: {sequence[4]} = 2^(5+1) + (3*{b_prev}-2) = 64 + {b}")

    print("\nStep 2: Calculating the next term in the sequence (n=6).")
    
    # Calculate next b
    b_prev = b
    final_b = 3 * b_prev - 2
    
    # Calculate next term
    power_of_2 = 2**(6+1)
    result = power_of_2 + final_b

    print("The final calculation is based on the discovered pattern:")
    print(f"{power_of_2} + {final_b} = {result}")

solve_sequence()