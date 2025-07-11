import math

def solve_message_length():
    """
    Calculates the maximum possible length of a secret message.
    """
    # Total number of symbols available in the Ching book.
    total_symbols = 9999

    # Determine the number of unique symbols that can be used for encoding.
    # There are 2 yinyang states (yin, yang) and 5 wuxing elements.
    yinyang_states = 2
    wuxing_elements = 5
    unique_encoding_symbols = yinyang_states * wuxing_elements

    # The number of unique characters the message can be written in.
    num_unique_characters = 120

    print(f"Step 1: Determine the parameters of the encoding problem.")
    print(f"The total number of symbols available in the Ching is {total_symbols}.")
    print(f"The number of unique encoding symbols is {yinyang_states} (yinyang) * {wuxing_elements} (wuxing) = {unique_encoding_symbols}.")
    print(f"The message is composed of a set of {num_unique_characters} unique characters.\n")

    # Determine the number of symbols from the Ching required to represent
    # one character from the Shu. Let this be L.
    # We need to find the smallest integer L such that unique_encoding_symbols^L >= num_unique_characters
    # i.e., 10^L >= 120
    # This is equivalent to L >= log10(120).
    symbols_per_character = math.ceil(math.log(num_unique_characters, unique_encoding_symbols))

    print(f"Step 2: Calculate the number of symbols needed to encode one character.")
    print(f"We need to find the smallest integer L such that {unique_encoding_symbols}^L is at least {num_unique_characters}.")
    print(f"For L=2, {unique_encoding_symbols}^2 = {unique_encoding_symbols**2}, which is less than {num_unique_characters}.")
    print(f"For L=3, {unique_encoding_symbols}^3 = {unique_encoding_symbols**3}, which is greater than {num_unique_characters}.")
    print(f"Therefore, we need {symbols_per_character} symbols from the Ching to represent one character of the message.\n")

    # Calculate the maximum possible length of the message.
    # This is the total number of symbols divided by the symbols needed per character.
    # We use integer division since we can't form a partial character.
    max_length = total_symbols // symbols_per_character

    print(f"Step 3: Calculate the maximum message length.")
    print(f"The final equation is the total symbols divided by the symbols per character:")
    print(f"{total_symbols} // {symbols_per_character} = {max_length}\n")
    print(f"The maximum length of the secret message is {max_length} characters.")

    return max_length

# Run the function to solve the problem
final_answer = solve_message_length()
print(f"<<<{final_answer}>>>")
