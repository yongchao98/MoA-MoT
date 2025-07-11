def solve_puzzle():
    """
    This function calculates the length of the secret word based on numbers
    extracted from the image and explains the reasoning to find the word.
    """
    # Numbers extracted from the image
    num_symbols_top = 22  # From the top border
    num_stars = 2         # Yellow stars
    num_cake_spots = 4    # Red spots on the cake

    # The equation to find the length of the secret word
    word_length = num_symbols_top / num_stars - num_cake_spots

    # Printing the equation step by step as requested
    print("The equation is derived from the elements in the image:")
    print(f"Number of symbols on top border: {num_symbols_top}")
    print(f"Number of stars: {num_stars}")
    print(f"Number of red spots on the cake: {num_cake_spots}")
    print("\nThe equation to find the word's length is:")
    print(f"{num_symbols_top} / {num_stars} - {num_cake_spots} = {int(word_length)}")
    
    print("\nThe result is 7. The character's clothing is the 'missing texture' pattern from Minecraft.")
    print("The 7-letter secret word is TEXTURE.")

solve_puzzle()