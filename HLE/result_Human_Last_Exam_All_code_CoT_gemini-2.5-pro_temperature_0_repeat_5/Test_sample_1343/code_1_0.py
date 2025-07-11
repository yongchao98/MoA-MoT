def solve_raspy_puzzle():
    """
    This function solves the puzzle by simulating the logic of the provided RASPy code.
    """

    def run_v_simulation(sop: str) -> str:
        """
        Simulates the logic of the RASPy function `v`.

        The core logic is a conditional check: if the first three characters of the
        input string contain a '7', a "pwned" message is returned. Otherwise,
        it performs addition on the numbers around the '+' sign.
        """
        # The condition is whether '7' is in the first three characters.
        has_seven_in_prefix = '7' in sop[:3]

        if has_seven_in_prefix:
            # This branch corresponds to the "get pwned" output.
            length = len(sop)
            base_msg = "get pwned"
            
            # The message is padded with '!' to match the input length.
            pwned_msg_list = list(base_msg + '!' * (length - len(base_msg)))
            
            # A special rule in the RASPy code changes the last character to '1'
            # if the input length is greater than 11.
            if length > 11:
                pwned_msg_list[length - 1] = '1'
            
            return "".join(pwned_msg_list)
        else:
            # This branch performs the addition.
            try:
                parts = sop.split('+')
                num1 = int(parts[0])
                num2 = int(parts[1])
                result = num1 + num2
                return str(result)
            except (ValueError, IndexError):
                # This is a fallback for malformed input.
                return "Invalid input for addition"

    # The two inputs for the puzzle
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    # --- Process Input 1 ---
    output1 = run_v_simulation(input1)

    # --- Process Input 2 ---
    # For the second input, we need to construct the full equation string.
    parts2 = input2.split('+')
    num1_str2 = parts2[0]
    num2_str2 = parts2[1]
    result2_str = run_v_simulation(input2)
    
    # The final output for the second part is the equation itself.
    output2_equation = f"{num1_str2} + {num2_str2} = {result2_str}"

    # Combine the results as per the "output1;output2" format.
    final_output = f"{output1};{output2_equation}"
    
    print(final_output)

solve_raspy_puzzle()