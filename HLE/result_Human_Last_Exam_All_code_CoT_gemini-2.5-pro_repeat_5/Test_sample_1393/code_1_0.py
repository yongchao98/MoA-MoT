import collections

def solve_puzzle():
    """
    This function decodes the Baudot-encoded answer choices and identifies the correct one.
    """
    
    # Baudot (ITA2) Letter Shift (LTRS) character map.
    baudot_ltrs_map = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10011': 'W', '11101': 'X', '10101': 'Y',
        '10001': 'Z', '00100': ' '
    }

    # The encoded answer choices provided in the puzzle.
    choices = {
        "A": "01111 00111 01100 00100 10111 00111",
        "B": "10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110",
        "C": "10101 00111 00001 00100 01011 00111",
        "D": "10101 00111 00100 01011 00111",
        "E": "01110 10100 00111 00001 01100 00100 01011 00111"
    }

    decoded_choices = {}
    print("--- Decoding Baudot Answer Choices ---")
    for label, code_string in choices.items():
        codes = code_string.split(' ')
        decoded_text = "".join([baudot_ltrs_map.get(code, '?') for code in codes])
        decoded_choices[label] = decoded_text
        print(f"{label}: {code_string} -> {decoded_text}")
    
    print("\n--- Analysis ---")
    print("The quote is attributed to the character Lü Bu. His boastful and dramatic personality is best captured by Chuan Opera (CHUAN JU), famous for its 'face-changing' technique.")
    print("Therefore, the correct answer is E.")

    # Per the instructions, showing the step-by-step 'equation' for the final answer.
    print("\n--- Detailed Decoding of Correct Answer (E) ---")
    correct_label = 'E'
    correct_codes = choices[correct_label].split(' ')
    
    equation_parts = []
    for code in correct_codes:
        char = baudot_ltrs_map.get(code, '?')
        equation_parts.append(f"{code} -> '{char}'")
    
    print(" + ".join(equation_parts))
    print(f"Final Result: {decoded_choices[correct_label]}")

solve_puzzle()