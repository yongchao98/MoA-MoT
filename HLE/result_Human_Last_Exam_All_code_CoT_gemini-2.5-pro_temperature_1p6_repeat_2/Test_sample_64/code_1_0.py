def find_pitch_accent_for_otouto():
    """
    This script finds and explains the Japanese pitch accent for the word 「弟」.
    """
    # Step 1: Define a dictionary to store pitch accent information.
    # This acts as our knowledge base for this query.
    # The notation [0] represents the Heiban pattern.
    word_database = {
        "弟": {
            "reading": "おとうと",
            "pattern_name": "Heiban (平板)",
            "notation": "[0]",
            "explanation": "Starts low, rises on the second mora, and stays high."
        }
    }

    target_word = "弟"
    
    # Step 2: Look up the word in our database.
    if target_word in word_database:
        info = word_database[target_word]
        
        print(f"Looking for the pitch accent of 「{target_word}」...")
        print(f"The standard pitch accent pattern is: {info['pattern_name']}")
        print(f"This pattern is also represented by the notation: {info['notation']}")
        
        # Step 3: Explain and visualize the pattern.
        # The word おとうと has 3 morae: お, とう, と
        morae = ["お (o)", "とう (tō)", "と (to)"]
        # Heiban [0] pattern is Low-High-High for a 3-mora word.
        pitch_visualization = "L-H-H  (L=Low, H=High)"
        
        print("\nPitch progression for おとうと:")
        print("Morae:   " + "   ".join(morae))
        print("Pitch:   " + pitch_visualization)
        
        print("\nBased on this information, the correct answer choice is A.")

    else:
        print(f"Pitch accent information for 「{target_word}」 not found.")

# Execute the function to find and display the answer.
find_pitch_accent_for_otouto()