def find_beatles_song():
    """
    This function identifies a Beatles song based on a specific chord progression.
    The user is asking for a song that features a chord jump from the tonic (I)
    to the minor fifth chord (vm) at the start.
    """

    # 1. Define the musical concepts.
    # The "tonic" is the first chord in a scale, the I chord.
    # The "minor fifth" is a minor chord built on the fifth note of the scale, the vm chord.
    # This progression is distinctive because the standard V chord is usually major.
    
    # 2. Analyze the song "Michelle".
    # The song "Michelle" is in the key of F Major.
    # Key: F Major
    # Tonic (I) chord: F Major
    # Fifth degree of the scale: C
    # Minor fifth (vm) chord: C minor

    # 3. Pinpoint the progression in the song.
    # The main vocal line begins "Michelle, ma belle...".
    # The chord progression for this line is:
    # "Michelle" -> F Major (The Tonic - I)
    # "ma belle" -> C minor (The Minor Fifth - vm)

    song_name = "Michelle"
    tonic_chord = "F Major"
    minor_fifth_chord = "C minor"

    print(f"The song by The Beatles that starts its main phrase with a jump from the tonic to the minor fifth chord is '{song_name}'.")
    print("\nHere is the musical analysis:")
    print(f"Key of the song: F Major")
    
    print("\nThe first chord of the verse ('Michelle...') is the tonic (I):")
    print(f"Chord 1: {tonic_chord}")
    
    print("\nThe second chord of the verse ('ma belle...') is the minor fifth (vm):")
    print(f"Chord 2: {minor_fifth_chord}")

find_beatles_song()