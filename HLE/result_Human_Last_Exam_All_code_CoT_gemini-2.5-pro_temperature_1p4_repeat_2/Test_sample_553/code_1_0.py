def solve_beatles_chord_puzzle():
    """
    This script identifies the Beatles song with a prominent jump from the
    tonic (I) chord to the minor fifth (v) chord. It breaks down the music
    theory behind the progression.
    """
    
    # 1. The song in question is "Michelle". It's in the key of F major.
    song_title = "Michelle"
    key = "F Major"
    scale_notes = ["F", "G", "A", "B♭", "C", "D", "E"]
    
    # 2. The "tonic" is the first degree of the scale. The chord built on it is the I chord.
    tonic_degree_numeral = "I"
    tonic_root_note = scale_notes[0]
    tonic_chord_name = f"{tonic_root_note} major"
    
    # 3. The "fifth" is the fifth degree of the scale. A "minor fifth chord"
    # is a minor chord built on this note, making it a "v" chord (lowercase).
    # This is a "borrowed" chord, as the natural chord on the fifth degree in a major key is major (V).
    fifth_degree_numeral = "v"
    fifth_root_note = scale_notes[4]
    minor_fifth_chord_name = f"{fifth_root_note} minor"
    
    # 4. Print the explanation and the final answer.
    print(f"The song by The Beatles that famously uses a chord jump from the tonic to the minor fifth is '{song_title}'.")
    print(f"\nThis progression occurs in the key of {key}.")
    
    print(f"\nThe tonic chord (represented by the Roman numeral '{tonic_degree_numeral}') is {tonic_chord_name}.")
    print(f"The minor fifth chord (represented by the Roman numeral '{fifth_degree_numeral}') is {minor_fifth_chord_name}.")
    
    print("\nWhile the song doesn't start with this exact jump, it's a defining feature of its bridge section ('Until I find a way...').")

    print("\nThe final 'equation' for the chord jump is:")
    print(f"{tonic_degree_numeral} -> {fifth_degree_numeral}")

# Execute the function to solve the puzzle.
solve_beatles_chord_puzzle()