def solve():
    """
    Solves the problem by simulating the logic of the provided RASPy code for the two given inputs.
    """
    
    def get_output_for_sop(sop: str) -> str:
        """
        This function determines the output for a single input string 'sop' by analyzing
        the trigger condition hidden in the RASPy function 'u' and then simulating
        the main function 'v'.
        """
        # The trigger in function 'u' depends on finding the digit '7' in three specific slices of the input.
        # q(sop) corresponds to the first 3 characters.
        s1 = sop[0:3]
        # r(sop) corresponds to characters at indices 3, 4, 5.
        s2 = sop[3:6]
        # p(sop) corresponds to the last 3 characters of the 15-char string (indices 12, 13, 14).
        s3 = sop[12:15]

        # The 'pwned' message is triggered if and only if '7' is present in all three slices.
        is_triggered = ('7' in s1) and ('7' in s2) and ('7' in s3)

        if is_triggered:
            # If the trigger is met, function 'v' outputs a hardcoded message.
            # The message is constructed to be "get pwned" followed by punctuation to match the input length,
            # with the very last character being '1'.
            return "get pwned!!!!!1"
        else:
            # If the trigger is not met, the program performs addition.
            try:
                parts = sop.split('+')
                num1 = int(parts[0])
                num2 = int(parts[1])
                result = num1 + num2
                return str(result)
            except (ValueError, IndexError):
                # Fallback for unexpected format, though not needed for the given inputs.
                return "Error: Invalid input format"

    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = get_output_for_sop(input1)
    output2 = get_output_for_sop(input2)

    # The problem asks for the two outputs to be joined by a semicolon.
    print(f"{output1};{output2}")

solve()