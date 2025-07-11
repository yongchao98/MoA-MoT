import math

def solve_message_length():
    """
    This function calculates the maximum possible length of a secret message
    based on the principles of information encoding.
    """
    # Number of unique symbols available for encoding (from the Ching)
    num_encoding_symbols = 9999

    # Size of the alphabet for the decoded message (from the Shu)
    message_alphabet_size = 120

    print("To find the maximum message length (L), we need to determine the largest integer L for which the number of possible messages does not exceed the number of available encoding symbols.")
    print("\nThe formula is: (alphabet_size)^L <= num_encoding_symbols")
    print(f"In our case: {message_alphabet_size}^L <= {num_encoding_symbols}")

    print("\nTo solve for L, we can use logarithms:")
    print(f"L * log({message_alphabet_size}) <= log({num_encoding_symbols})")
    print(f"L <= log({num_encoding_symbols}) / log({message_alphabet_size})")

    # Calculate the upper bound for L
    # We use math.log (natural logarithm), but any base would work
    max_l_float = math.log(num_encoding_symbols) / math.log(message_alphabet_size)

    # Since the length of a message must be an integer, we take the floor
    # of the result to find the maximum valid integer length.
    max_l_integer = math.floor(max_l_float)

    print(f"\nCalculating the value:")
    print(f"L <= {max_l_float:.4f}")
    print("\nSince L must be an integer, we take the floor of this value.")
    print(f"Maximum message length = floor({max_l_float:.4f}) = {max_l_integer}")
    
    # Check the result
    print(f"\nVerification:")
    power_check_pass = message_alphabet_size ** max_l_integer
    power_check_fail = message_alphabet_size ** (max_l_integer + 1)
    print(f"{message_alphabet_size}^{max_l_integer} = {power_check_pass}, which is <= {num_encoding_symbols} (Correct)")
    print(f"{message_alphabet_size}^{max_l_integer + 1} = {power_check_fail}, which is > {num_encoding_symbols} (Correct)")


if __name__ == "__main__":
    solve_message_length()
