def find_beatles_song():
    """
    Identifies the Beatles song based on a specific harmonic progression
    and explains the music theory behind it.
    """
    song_title = "Michelle"
    key = "F minor"
    first_chord = "C major"
    second_chord = "C minor"
    
    # In music theory, Roman numerals are used to denote chords within a key.
    # 'V' represents the major chord built on the fifth degree of the scale (the dominant).
    # 'v' represents the minor chord built on the fifth degree of the scale.
    first_chord_function = "V"
    second_chord_function = "v"

    print(f"The song by The Beatles that features this distinctive chord movement is '{song_title}'.")
    print("\nHere's the music theory breakdown:")
    print(f"1. The song's verse section is in the key of {key}.")
    print(f"2. The famous opening vocal line ('Michelle, ma belle...') begins over a {first_chord} chord, which then 'jumps' to a {second_chord} chord.")
    print(f"3. In the key of {key}, the C chord is the fifth (the dominant).")
    print("\nTherefore, the progression is from the major fifth to the minor fifth.")
    print("The final 'equation' for this harmonic progression is:")
    print(f"{first_chord_function} -> {second_chord_function}")

if __name__ == "__main__":
    find_beatles_song()