def solve_raspy_puzzle():
    """
    Solves the RASPy puzzle by implementing the deduced high-level logic.
    """
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    def u_trigger_logic(sop: str) -> int:
        """
        This function simulates the logic of the RASPy function 'u'.
        It returns 1 if the '7' character is present in three specific
        3-character chunks of the input string, and 0 otherwise.
        """
        # Chunk 1: First 3 characters (from function q)
        chunk1 = sop[0:3]
        # Chunk 2: Characters at indices 3, 4, 5 (from function r)
        chunk2 = sop[3:6]
        # Chunk 3: Last 3 characters (from function p)
        chunk3 = sop[-3:]

        has_seven_in_chunk1 = '7' in chunk1
        has_seven_in_chunk2 = '7' in chunk2
        has_seven_in_chunk3 = '7' in chunk3

        if has_seven_in_chunk1 and has_seven_in_chunk2 and has_seven_in_chunk3:
            return 1
        else:
            return 0

    def v_logic(sop: str) -> str:
        """
        This function simulates the logic of the main RASPy function 'v'.
        It branches based on the result of the 'u' function logic.
        """
        trigger = u_trigger_logic(sop)

        if trigger == 1:
            # If the trigger is 1, return the "pwned" message.
            # The "aethetics" part in the original code changes the last character.
            return "get pwned!!!1"
        else:
            # If the trigger is 0, perform the addition.
            try:
                parts = sop.split('+')
                num1 = int(parts[0])
                num2 = int(parts[1])
                result = num1 + num2
                # The prompt asks to omit leading zeroes, which str() does by default.
                return str(result)
            except (ValueError, IndexError):
                return "Invalid input format"

    # Calculate the output for both inputs
    output1 = v_logic(input1)
    output2 = v_logic(input2)

    # Print the final result in the specified format "output1;output2"
    print(f"{output1};{output2}")

solve_raspy_puzzle()