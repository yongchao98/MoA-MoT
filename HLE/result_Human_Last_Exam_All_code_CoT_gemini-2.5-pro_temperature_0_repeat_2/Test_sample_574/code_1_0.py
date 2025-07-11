def solve_ship_name_puzzle():
    """
    This script solves a word puzzle to find a ship name from the Culture series.
    """
    
    # Step 1: The source word for "sleeveless garments..." is "BOLEROS".
    # From "BOLEROS", we can form the word "SOBER".
    word1 = "SOBER"
    
    # Step 2: The source word for "experienced and trusted individuals..." is "COUNSELORS".
    # From "COUNSELORS", we can form the word "COUNSEL".
    word2 = "COUNSEL"
    
    # Step 3: Combine the two words to get the ship name.
    ship_name = f"{word1} {word2}"
    
    # Step 4: Print the final equation, showing each word used.
    print("The puzzle is solved as follows:")
    print(f"The first word is '{word1}'.")
    print(f"The second word is '{word2}'.")
    print("The combined ship name is derived from the equation:")
    print(f"{word1} + {word2} = {ship_name}")

solve_ship_name_puzzle()