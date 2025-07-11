def solve_puzzle():
    """
    This script decodes a Baudot-encoded message to find the origin of a famous quote.
    """

    # International Telegraph Alphabet No. 2 (ITA2) LTRS mapping
    ita2_map = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10011': 'W', '11101': 'X', '10101': 'Y',
        '10001': 'Z', '00100': ' '
    }

    # The encoded message for choice E
    encoded_choice_e = ['01110', '10100', '00111', '00011', '01100', '00100', '01011', '00111']
    
    decoded_message = []
    print("Decoding Choice E:")
    
    # Decode each part of the message and print the process
    for code in encoded_choice_e:
        char = ita2_map.get(code, '?')
        print(f"{code} -> {char}")
        decoded_message.append(char)
        
    final_answer = "".join(decoded_message)
    print("\nFinal Decoded Answer:")
    print(final_answer)

solve_puzzle()
<<<E>>>