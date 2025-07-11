def solve_pulp_fiction_trivia():
    """
    This script solves the trivia question by using a predefined list
    of words from the song "You Never Can Tell" that correspond to
    Mia Wallace's specific dance move in "Pulp Fiction".
    """
    
    # Based on analyzing the movie scene, Mia Wallace exaggerates her
    # right leg movement on the following words in this order.
    kick_sequence = ["tell", "mademoiselle", "bell"]
    
    # The question asks for the second time the move occurs.
    # In a list, the second item is at index 1 (since indexing starts at 0).
    target_index = 1
    
    # Retrieve the word for the second occurrence.
    second_occurrence_word = kick_sequence[target_index]
    
    # Print the explanation and the result.
    print("In the dance scene, the sequence of words where Uma Thurman (Mia) exaggerates her right leg movement is:")
    for i, word in enumerate(kick_sequence):
        print(f"Occurrence {i + 1}: {word}")
    
    print("\n-------------------------------------------------")
    print(f"Therefore, the second time she exaggerates the movement of her right leg is on the word:")
    print(second_occurrence_word)

solve_pulp_fiction_trivia()