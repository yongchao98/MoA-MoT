def solve_raspy_puzzle():
    """
    This function simulates the logic of the provided RASPy code to determine the output
    for the two given inputs.
    """

    def calculate_final_output(input_string):
        """
        Calculates the output of the v function for a given input string.
        Based on the analysis, the u() function always triggers the "get pwned" branch.
        The final string's content depends only on the input string's length.
        """
        length = len(input_string)

        # The base message is "get pwned" converted to ASCII ordinal values.
        base_msg_ord = [103, 101, 116, 32, 112, 119, 110, 101, 100]

        # The message is padded with '!' (ASCII 33) to match the input length.
        num_padding_chars = length - len(base_msg_ord)
        
        if num_padding_chars >= 0:
            output_ord_list = base_msg_ord + [33] * num_padding_chars
        else:
            # In case the input string is shorter than "get pwned"
            output_ord_list = base_msg_ord[:length]

        # A rule modifies the last character to '1' (ASCII 49) if the index is > 10.
        # This corresponds to the string length being > 11.
        last_index = length - 1
        if last_index > 10:
            output_ord_list[last_index] = 49

        # Convert the list of ASCII values back to a character string.
        final_string = "".join([chr(c) for c in output_ord_list])
        return final_string

    # First input string
    input1 = "734107+4295754"
    output1 = calculate_final_output(input1)

    # Second input string
    input2 = "5429141+142196"
    output2 = calculate_final_output(input2)

    # Print the results in the specified "output1;output2" format.
    print(f"{output1};{output2}")

solve_raspy_puzzle()