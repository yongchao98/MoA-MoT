import struct

def solve_problem():
    """
    Analyzes the C program and generates the final answer in the 'ab' format.
    """

    # Part a: Is the program correct on the two given examples?
    # Simulate the C code's logic to verify correctness.

    def simulate_word_abbreviation(input_word):
        """
        Simulates the C program's logic for a given word.
        """
        # unsigned long long s is treated as an 8-byte buffer.
        s_buffer = bytearray(8)
        word_len = 0

        # Simulate the 'add' function
        for char in input_word:
            if word_len < 7:
                s_buffer[word_len] = ord(char)
            else:
                s_buffer[7] = ord(char)
            word_len += 1

        # Simulate the main conditional and print functions
        if word_len > 8:
            # printl logic
            first_char = chr(s_buffer[0])
            middle_num = word_len - 2
            last_char = chr(s_buffer[7])
            return f"{first_char}{middle_num}{last_char}", s_buffer
        else:
            # prints logic
            return s_buffer[:word_len].decode('ascii'), s_buffer

    # Check example 1: "localization"
    output_loc, s_bytes_loc = simulate_word_abbreviation("localization")
    is_correct_loc = (output_loc == "l10n")

    # Check example 2: "internationalization"
    output_int, _ = simulate_word_abbreviation("internationalization")
    is_correct_int = (output_int == "i18n")

    answer_a = 'Y' if is_correct_loc and is_correct_int else 'N'

    # Part b: Is it correct for every input? Otherwise answer value of s.
    # The program's logic is sound on a typical little-endian machine.
    # Therefore, we provide the value of 's' for the input "localization".
    
    # The bytes stored in 's' for "localization" are in s_bytes_loc
    # s_bytes_loc contains the ASCII values of 'l','o','c','a','l','i','z','n'
    # We convert these 8 bytes to an unsigned 64-bit integer (ull),
    # assuming little-endian byte order as is standard.
    s_value = struct.unpack('<Q', s_bytes_loc)[0]
    
    # Format the value as a hexadecimal string for the answer.
    answer_b = hex(s_value)

    # Combine answers into the final 'ab' format and print.
    final_answer = f"{answer_a}{answer_b}"
    print(final_answer)

solve_problem()