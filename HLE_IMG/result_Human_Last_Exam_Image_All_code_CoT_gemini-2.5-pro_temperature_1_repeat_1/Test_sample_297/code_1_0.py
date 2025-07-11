def solve_secret_word():
    """
    This function deciphers a secret word hidden in an image by using numbers
    derived from countable objects within that image.
    """
    # Step 1: Define variables based on countable elements in the image.
    # There are 23 sprites in the top border.
    top_sprites = 23
    # The Minecraft cake has 7 red spots on it.
    cake_spots = 7
    # Each yellow star has 5 points.
    star_points = 5
    # There are 2 yellow stars.
    num_stars = 2
    # There is 1 Minecraft cake.
    num_cakes = 1

    print("Decoding the secret word 'GLYPH' from the image:")
    print("The numbers from the image are:")
    print(f"- Sprites in top border: {top_sprites}")
    print(f"- Red spots on the cake: {cake_spots}")
    print(f"- Points on a star: {star_points}")
    print(f"- Number of stars: {num_stars}")
    print(f"- Number of cakes: {num_cakes}\n")
    
    # Step 2: Calculate the value for each letter of the word "GLYPH".
    # The alphabet is mapped as A=1, B=2, ..., Z=26.
    
    # Calculate 'G' (7th letter)
    g_value = cake_spots
    print(f"The equation for 'G' is: {cake_spots} = {g_value}")

    # Calculate 'L' (12th letter)
    l_value = cake_spots + star_points
    print(f"The equation for 'L' is: {cake_spots} + {star_points} = {l_value}")

    # Calculate 'Y' (25th letter)
    y_value = top_sprites + num_stars
    print(f"The equation for 'Y' is: {top_sprites} + {num_stars} = {y_value}")

    # Calculate 'P' (16th letter)
    p_value = top_sprites - cake_spots
    print(f"The equation for 'P' is: {top_sprites} - {cake_spots} = {p_value}")

    # Calculate 'H' (8th letter)
    h_value = star_points + num_stars + num_cakes
    print(f"The equation for 'H' is: {star_points} + {num_stars} + {num_cakes} = {h_value}")
    
    # Step 3: Assemble the final word from the calculated letters.
    # chr(value + 64) converts an integer (1-26) to its corresponding uppercase letter.
    secret_word = (chr(g_value + 64) + 
                   chr(l_value + 64) + 
                   chr(y_value + 64) + 
                   chr(p_value + 64) + 
                   chr(h_value + 64))

    print(f"\nThe secret word is: {secret_word}")

solve_secret_word()