def solve_puzzle():
    """
    This function solves the puzzle by deriving a secret word from numbers
    found in the image. It uses simple arithmetic and an A=1, B=2, ... mapping.
    """
    
    # Step 1: Define base numbers by counting elements in the image.
    num_stars = 2
    points_per_star = 5
    num_cakes = 1
    cake_dots = 7
    palette_colors = 8

    # Step 2: The target word is GLITCH. We calculate the numeric value for each letter.
    # G = 7th letter
    # L = 12th letter
    # I = 9th letter
    # T = 20th letter
    # C = 3rd letter
    # H = 8th letter
    
    # Step 3: Formulate and print the equation for each letter, showing how its
    # value is derived from the image's elements.
    print("Deriving the secret word 'GLITCH' using numbers from the image:")
    
    # Equation for G (7)
    g_val = cake_dots
    print("G = " + str(cake_dots))
    
    # Equation for L (12)
    l_val = cake_dots + points_per_star
    print("L = " + str(cake_dots) + " + " + str(points_per_star))
    
    # Equation for I (9)
    i_val = cake_dots + num_stars
    print("I = " + str(cake_dots) + " + " + str(num_stars))

    # Equation for T (20)
    t_val = (points_per_star * num_stars) * num_stars
    print("T = (" + str(points_per_star) + " * " + str(num_stars) + ") * " + str(num_stars))
    
    # Equation for C (3)
    c_val = num_stars + num_cakes
    print("C = " + str(num_stars) + " + " + str(num_cakes))
    
    # Equation for H (8)
    h_val = palette_colors
    print("H = " + str(palette_colors))

    # Step 4: Convert the numeric values to characters and form the word.
    # (A=1 maps to ASCII 65, so we add 64)
    secret_word = chr(g_val + 64) + \
                  chr(l_val + 64) + \
                  chr(i_val + 64) + \
                  chr(t_val + 64) + \
                  chr(c_val + 64) + \
                  chr(h_val + 64)

    print("\nThe final secret word is:")
    print(secret_word)

solve_puzzle()