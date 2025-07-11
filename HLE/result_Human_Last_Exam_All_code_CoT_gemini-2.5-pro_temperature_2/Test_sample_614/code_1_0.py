def solve_riddle():
    """
    This function solves the number riddle based on the provided clues
    and prints the final ordered sequence.
    """

    # The sequence is derived from the step-by-step logical deductions from the riddle's clues:
    # 1. Number 5 is last. (Sequence: _, _, _, _, 5)
    # 2. Number 3 protects the last. (Sequence: _, _, _, 3, 5)
    # 3. Number 1 is first. (Sequence: 1, _, _, 3, 5)
    # 4. Number 4 goes behind number 2, filling the remaining spots. (Sequence: 1, 2, 4, 3, 5)
    final_sequence = [1, 2, 4, 3, 5]

    # The prompt requests to output each number in the "final equation".
    # We will present the solution as a final, ordered sequence.
    num1, num2, num3, num4, num5 = final_sequence
    
    print(f"The final sequence is: {num1}, {num2}, {num3}, {num4}, {num5}")

solve_riddle()