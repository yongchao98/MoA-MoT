import string

def solve_nolan_cipher():
    """
    Solves a two-layer cipher (Substitution + Caesar) to find a quote
    from a Christopher Nolan movie.
    """
    # 1. SETUP: Define the clues from the prompt.
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    cipher_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    plain_alphabet = string.ascii_uppercase
    shift_key = 5

    # 2. DECODING LOGIC:
    # Hypothesis: The original message was Caesar-shifted (+5), then substituted.
    # To decode, we must reverse the substitution, then reverse the Caesar shift (-5).

    # Create a mapping to reverse the substitution cipher.
    # Maps from the custom cipher alphabet back to the plain A-Z alphabet.
    # e.g., 'C' -> 'A', 'H' -> 'B', 'R' -> 'C'
    reverse_sub_map = {cipher_char: plain_char for cipher_char, plain_char in zip(cipher_alphabet, plain_alphabet)}

    decoded_message = ""
    print("Decoding the message step-by-step:")
    print("-" * 35)

    for char in encoded_message:
        if char == ' ':
            decoded_message += ' '
            print("\n") # Add a line break for the space between words
            continue

        # Step A: Reverse the Substitution
        intermediate_char = reverse_sub_map[char]
        intermediate_ord = ord(intermediate_char)

        # Step B: Reverse the Caesar Shift
        # New index = (current index - shift amount + 26) % 26
        final_ord = ((intermediate_ord - ord('A') - shift_key) + 26) % 26 + ord('A')
        final_char = chr(final_ord)

        # Build the final decoded message string
        decoded_message += final_char

        # Print the equation for the current character
        print(f"'{char}' -> unsub to '{intermediate_char}' -> shift by -{shift_key} -> '{final_char}'")

    print("-" * 35)
    print(f"\nDecoded Message: {decoded_message}")
    print("\nThis famous quote is said by the character Eames in the movie Inception.")


solve_nolan_cipher()
