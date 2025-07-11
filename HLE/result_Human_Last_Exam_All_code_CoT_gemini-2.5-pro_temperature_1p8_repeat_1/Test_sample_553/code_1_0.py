def find_beatles_song():
    """
    This function identifies and explains the Beatles song
    that begins with a distinctive harmonic jump from the
    tonic chord to a minor chord on the fifth degree.
    """
    song_title = "Michelle"
    key = "F Major"
    
    # In music theory, chords are numbered with Roman numerals.
    # 'I' represents the tonic chord.
    # 'v' represents a minor chord built on the fifth degree of the scale.
    tonic_chord_name = "F Major"
    tonic_chord_numeral = "I"
    
    minor_fifth_chord_name = "C minor"
    minor_fifth_chord_numeral = "v"

    print(f"The song by The Beatles that fits the description is: '{song_title}'")
    print("\nHere is the musical explanation:")
    print(f"The song is in the key of {key}.")
    print(f"The opening vocal harmony starts on the tonic chord, which is '{tonic_chord_name}'.")
    print(f"It then immediately moves to a '{minor_fifth_chord_name}' chord, which is a minor chord built on the fifth degree of the scale.")
    
    print("\nThe harmonic 'equation' for this opening is:")
    print(f"Chord 1 (Tonic): {tonic_chord_numeral} = {tonic_chord_name}")
    print(f"Chord 2 (Minor Fifth): {minor_fifth_chord_numeral} = {minor_fifth_chord_name}")

if __name__ == "__main__":
    find_beatles_song()