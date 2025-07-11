import sys
import io

# Capture original stdout to restore it later if needed, though not strictly necessary here.
original_stdout = sys.stdout
# Create a string buffer to hold the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# Plan:
# 1. Define the Baudot to Character mapping (ITA2 LTRS shift).
# 2. Define the encrypted answer choices.
# 3. Decode each choice and print the process.
# 4. Identify the correct choice based on the context of the puzzle.
#    - The quote is from Zhang Fei, a general of Shu-Han (in modern Sichuan).
#    - The origin is likely a cultural form from that region, i.e., Sichuan Opera ("Chuan Ju").
# 5. Print the final answer with a clear explanation and the requested code breakdown.

def solve_baudot_puzzle():
    """
    Decodes the answer choices and identifies the correct origin of the quote.
    """
    # ITA2 Baudot Code Mapping (Letters Shift)
    baudot_map = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10011': 'W', '11101': 'X', '10101': 'Y',
        '10001': 'Z', '00100': ' '
    }

    # Answer choices encoded in Baudot code
    choices = {
        "A": "01111 00111 01100 00100 10111 00111",
        "B": "10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110",
        "C": "10101 00111 00001 00100 01011 00111",
        "D": "10101 00111 00100 01011 00111",
        "E": "01110 10100 00111 00011 01100 00100 01011 00111"
    }

    correct_choice_label = "E"
    correct_decoded_text = "CHUAN JU"

    print("Step 1: The decoded Morse code asks for the origin of a quote by Zhang Fei from 'Romance of the Three Kingdoms'.")
    print("Step 2: Zhang Fei is strongly associated with the Shu-Han kingdom, located in modern-day Sichuan.")
    print("Step 3: The answer choices are decoded from Baudot code below to find a cultural origin related to Sichuan.")
    print("-" * 30)

    for option, encoded_str in choices.items():
        codes = encoded_str.split(' ')
        decoded_word = "".join([baudot_map.get(code, "?") for code in codes])
        print(f"Decoding Option {option}: {decoded_word}")
    
    print("-" * 30)
    print("Explanation:")
    print(f"Option E decodes to 'CHUAN JU', which means Sichuan Opera.")
    print("This is the correct answer because Sichuan Opera ('Chuan Ju') is a famous cultural art form from the region historically associated with Zhang Fei, and it frequently features his stories.")
    
    print("\nFinal Answer Breakdown:")
    print("The codes for the correct answer 'CHUAN JU' are:")
    
    codes = choices[correct_choice_label].split(' ')
    letters = list(correct_decoded_text)
    
    output_parts = []
    for i in range(len(codes)):
        char = letters[i]
        code = codes[i]
        if char == ' ':
             output_parts.append(f"SPACE({code})")
        else:
             output_parts.append(f"{char}({code})")
    
    # "output each number in the final equation!" -> This is interpreted as showing the full code mapping.
    print(" + ".join(output_parts))

solve_baudot_puzzle()

# Get the content from the buffer
final_output = output_buffer.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the captured output to the actual console
print(final_output)
print("<<<E>>>")