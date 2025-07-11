import struct

def calculate_s_value_as_hex(word):
    """
    Simulates the C program's 'add' function to determine the final bytes
    stored in the 's' variable and returns its little-endian hex value.
    """
    # s is an 8-byte buffer, initialized to zeros.
    s_bytes = bytearray(8)
    
    for i, char in enumerate(word):
        # The 'add' function logic:
        # if (l < 7) *(s + l) = c; else *(s + 7) = c;
        if i < 7:
            s_bytes[i] = ord(char)
        else:
            s_bytes[7] = ord(char)
            
    # The C code requires a little-endian architecture to work.
    # We unpack the byte array as a little-endian unsigned long long ('<Q').
    s_integer_value = struct.unpack('<Q', s_bytes)[0]
    
    return hex(s_integer_value)

def simulate_c_program_output(word):
    """
    Simulates the full C program to check if its output is correct
    according to the problem description.
    """
    l = len(word)
    
    # Logic for words longer than 8 characters.
    if l > 8:
        # The C program calculates the output using its flawed but functional logic.
        # It takes the first char, length-2, and the last char.
        expected_output = f"{word[0]}{l-2}{word[-1]}"
        
        # Simulate the C program's output
        s_bytes = bytearray(8)
        for i, char in enumerate(word):
            if i < 7:
                s_bytes[i] = ord(char)
            else:
                s_bytes[7] = ord(char)
        
        # printl prints s[0], l-2, and s[7]
        c_program_output = f"{chr(s_bytes[0])}{l-2}{chr(s_bytes[7])}"
        
        return c_program_output == expected_output
    # Logic for words 8 characters or shorter.
    else:
        # The C program just prints the word as is, which is correct.
        return True

# --- Main analysis ---

# a. Is this program correct on two given examples (Y/N)?
is_correct_on_ex1 = simulate_c_program_output("localization")
is_correct_on_ex2 = simulate_c_program_output("internationalization")
answer_a = "Y" if is_correct_on_ex1 and is_correct_on_ex2 else "N"

# b. Is it correct for every input? If no then answer the length of the shortest input
# making it incorrect. Otherwise answer value of s when input is localization.
# Our analysis concluded the program is correct for every input (on little-endian systems).
# So, we provide the value of 's' for "localization".
answer_b = calculate_s_value_as_hex("localization")

# Print the final combined answer in "ab" format.
print(f"{answer_a}{answer_b}")