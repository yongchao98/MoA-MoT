import string

def solve_puzzle():
    """
    This function solves the puzzle by calculating the indices of the secret word's
    letters using numbers found in the image.
    """
    # Step 1: Define the numbers found in the image.
    num_stars = 2
    num_cake_spots = 7
    num_border_sprites = 16

    print("Decoding the secret word using numbers from the image...")
    print(f"Number of stars: {num_stars}")
    print(f"Number of red spots on the cake: {num_cake_spots}")
    print(f"Number of sprites in the top border: {num_border_sprites}\n")

    # Step 2: Define the alphabet for mapping indices to letters (1-based index).
    alphabet = string.ascii_uppercase

    # Step 3: Calculate the index for each letter of the secret word.
    
    # Calculate index for 'C'
    # Equation: 7 - (16 / 2 / 2) = 3
    c_index = num_cake_spots - (num_border_sprites // 2 // 2)
    print(f"First letter's index: {num_cake_spots} - ({num_border_sprites} / 2 / 2) = {c_index}")

    # Calculate index for 'A'
    # Equation: (16 / 2) - 7 = 1
    a_index = (num_border_sprites // 2) - num_cake_spots
    print(f"Second letter's index: ({num_border_sprites} / 2) - {num_cake_spots} = {a_index}")

    # Calculate index for 'K'
    # Equation: 16 - (7 - 2) = 11
    k_index = num_border_sprites - (num_cake_spots - num_stars)
    print(f"Third letter's index: {num_border_sprites} - ({num_cake_spots} - {num_stars}) = {k_index}")
    
    # Calculate index for 'E'
    # Equation: 7 - 2 = 5
    e_index = num_cake_spots - num_stars
    print(f"Fourth letter's index: {num_cake_spots} - {num_stars} = {e_index}\n")

    # Step 4: Retrieve the letters and construct the final word.
    try:
        char_c = alphabet[c_index - 1]
        char_a = alphabet[a_index - 1]
        char_k = alphabet[k_index - 1]
        char_e = alphabet[e_index - 1]
        secret_word = char_c + char_a + char_k + char_e
        
        print(f"The calculated indices ({c_index}, {a_index}, {k_index}, {e_index}) map to the letters '{char_c}', '{char_a}', '{char_k}', '{char_e}'.")
        print(f"\nThe secret word is: {secret_word}")

    except IndexError:
        print("An error occurred while mapping indices to letters.")

if __name__ == '__main__':
    solve_puzzle()
