def find_final_note():
    """
    This function explains the logic for finding the final sung note
    in the provided arrangement of "Happy Birthday to You".
    """

    # Step 1: Determine the chord for the final "you".
    print("Step 1: Mapping the chords to the song's structure.")
    print("The problem states a consistent pattern:")
    print("- 'birthday' is accompanied by [X]m7 followed by a 7(9) chord.")
    print("- 'you' is accompanied by a chord that is a half-step down from [X]m7.")
    print("\nSong Structure Analysis:")
    print("1. Phrase 1 ('...to you'): Starts with Cm7, so 'you' is on Bm7.")
    print("2. Phrase 2 ('...to you'): Starts with Am7, so 'you' is on Abm7.")
    print("3. Phrase 3 ('...dear Name'): Starts with Ebm7, so 'Name' is on Bm7 (as per the provided chord list).")
    print("4. Phrase 4 (Final '...to you'): The progression restarts with Cm7, so the final 'you' is on Bm7.")
    final_chord = "Bm7"
    print(f"\nResult of Step 1: The chord played during the final 'you' is {final_chord}.")
    print("-" * 30)

    # Step 2: Determine the key of the song.
    print("Step 2: Establishing the song's key (tonal center).")
    starting_chord = "Cm7"
    print(f"The chord progression begins and cyclically ends with {starting_chord}.")
    print("This establishes the key of the song as C.")
    song_key = "C"
    print(f"\nResult of Step 2: The key of the song is {song_key}.")
    print("-" * 30)

    # Step 3: Identify the final melody note.
    print("Step 3: Identifying the final note of the melody.")
    print("The traditional melody of 'Happy Birthday to You' always resolves to the tonic note of the key on the final word.")
    final_note = song_key
    print(f"The tonic is the root of the key. Since the key is C, the final sung note must be C.")
    print("-" * 30)
    
    # Step 4: Final Conclusion
    # The requested format is to show the equation.
    print("Final Conclusion:")
    print(f"The concluding sung note of the song is the tonic of its key.")
    print(f"Final Sung Note = {final_note}")

find_final_note()