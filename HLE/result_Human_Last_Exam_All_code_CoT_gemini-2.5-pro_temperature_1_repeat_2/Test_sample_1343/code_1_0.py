def solve():
    """
    This function simulates the logic of the provided RASPy code for the two given inputs.
    """
    
    inputs = ["734107+4295754", "5429141+142196"]
    outputs = []

    for sop in inputs:
        # Step 1: Analyze the trigger condition from function u(sop).
        # The trigger requires the digit '7' to be present in sop[0:3], sop[3:6], and sop[-3:].
        part1 = sop[0:3]
        part2 = sop[3:6]
        part3 = sop[-3:]
        
        trigger = ('7' in part1) and ('7' in part2) and ('7' in part3)

        if trigger:
            # Step 2: If triggered, produce the "get pwned" message.
            length = len(sop)
            # This corresponds to the line: ee = [103, 101, 116, 32, 112, 119, 110, 101, 100] + [33] * 89
            base_msg = "get pwned"
            # Pad with '!' to match the full length of the input.
            result_str = base_msg.ljust(length, '!')
            
            # This corresponds to the "aesthetics" line, which changes the last character to '1'
            # if the length is appropriate (which it is for the first input).
            # (dd == 1) & (indices + 1 == length) & (indices > 10)
            # For len=14, this is true at index 13.
            if length > 11:
                 result_list = list(result_str)
                 result_list[length - 1] = '1'
                 result_str = "".join(result_list)
            
            outputs.append(result_str)
        else:
            # Step 3: If not triggered, calculate the sum as specified.
            # This simulates the path that calculates aa, bb, and cc.
            
            # The prompt requires printing each number in the equation.
            # We will first calculate the sum and then format the output.
            parts = sop.split('+')
            num1_str = parts[0]
            num2_str = parts[1]

            # The calculation in RASPy is equivalent to standard integer addition.
            result = int(num1_str) + int(num2_str)
            
            # Format the output as requested by the problem.
            equation = f"{num1_str} + {num2_str} = {result}"
            # The final output requested for the prompt is just the result number.
            # We will store the number, and print the equation for clarity as requested.
            print(f"Calculating sum for '{sop}': {equation}")
            outputs.append(str(result))

    # Step 4: Print the final results in the specified format "output1;output2".
    print("\nFinal combined output:")
    final_output_string = ";".join(outputs)
    print(final_output_string)

solve()
<<<get pwned!!!1;5571337>>>