def solve():
    """
    This function implements the logic derived from analyzing the RASPy code
    to solve the problem for both inputs.
    """

    def check_trigger_condition(input_string):
        """
        This function simulates the logic of the RASPy function 'u'.
        It checks if the digit '7' is present in three specific chunks
        of the input string.
        - Chunk 1: First 3 characters.
        - Chunk 2: Characters at indices 3 through 5.
        - Chunk 3: Last 3 characters.
        The trigger condition is met only if all three chunks contain a '7'.
        """
        # Define the chunks based on the RASPy functions q, r, and p.
        chunk1 = input_string[0:3]
        chunk2 = input_string[3:6]
        chunk3 = input_string[-3:]

        # Check for the presence of '7' in each chunk.
        cond1 = '7' in chunk1
        cond2 = '7' in chunk2
        cond3 = '7' in chunk3

        return cond1 and cond2 and cond3

    # --- Process Input 1 ---
    input1 = "734107+4295754"
    if check_trigger_condition(input1):
        # The length of the output string is the same as the input string (15).
        # "get pwned" (9 chars) + "!!!!!!" (6 chars) -> "get pwned!!!!!!"
        # The last char '!' is replaced by '1'.
        output1 = "get pwned!!!!!1"
    else:
        # This branch is not taken for input 1.
        parts = input1.split('+')
        num1 = int(parts[0])
        num2 = int(parts[1])
        output1 = str(num1 + num2)

    # --- Process Input 2 ---
    input2 = "5429141+142196"
    if check_trigger_condition(input2):
        # This branch is not taken for input 2.
        output2 = "get pwned... (error in logic)"
    else:
        parts = input2.split('+')
        num1 = int(parts[0])
        num2 = int(parts[1])
        result = num1 + num2
        # The problem asks to omit leading zeros, so we just use the standard string representation.
        output2 = str(result)
        # Note: The problem asked to "output each number in the final equation".
        # This is how the equation for the second input looks:
        # print(f"{num1} + {num2} = {result}")

    # --- Combine and print the final result in the specified format ---
    print(f"{output1};{output2}")

solve()
<<<get pwned!!!!!1;5571337>>>