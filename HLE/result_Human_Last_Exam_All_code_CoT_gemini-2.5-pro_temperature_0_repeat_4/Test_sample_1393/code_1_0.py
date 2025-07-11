def solve_baudot_puzzle():
    """
    This function decodes a Baudot-encoded message to find the origin of a quote.
    The quote is from the story of Zhao Yun (Three Kingdoms period), which has strong
    ties to the Sichuan region. We will decode option E, which corresponds to
    "CHUAN JU" (Sichuan Opera).
    """

    # Baudot (ITA2) LTRS shift character mapping
    ITA2_LTRS = {
        "00011": "A", "11001": "B", "01110": "C", "01001": "D", "00001": "E",
        "01101": "F", "11010": "G", "10100": "H", "00110": "I", "01011": "J",
        "01111": "K", "10010": "L", "11100": "M", "01100": "N", "11000": "O",
        "10110": "P", "11101": "Q", "01010": "R", "00101": "S", "10000": "T",
        "00111": "U", "11110": "V", "10011": "W", "11101": "X", # In ITA2, X can be the same as Q
        "10101": "Y", "10001": "Z",
        "00100": " "  # Space
    }

    # The Baudot code for option E
    option_e_code = "01110 10100 00111 00011 01100 00100 01011 00111"
    
    codes = option_e_code.split()
    decoded_chars = []

    print("Decoding the Baudot code for the correct option (E):")
    
    # This loop decodes each 5-bit number and prints the mapping
    for code in codes:
        char = ITA2_LTRS.get(code, "?")
        decoded_chars.append(char)
        # The problem asks to "output each number in the final equation"
        # We interpret this as showing the code-to-character mapping.
        print(f"{code} -> {char}")

    final_text = "".join(decoded_chars)
    print(f"\nThe fully decoded text is: {final_text}")

solve_baudot_puzzle()
<<<E>>>