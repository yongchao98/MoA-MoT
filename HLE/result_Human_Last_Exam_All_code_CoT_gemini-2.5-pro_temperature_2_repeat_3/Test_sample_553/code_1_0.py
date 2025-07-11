def find_beatles_song():
    """
    Identifies the Beatles song based on a description of its opening chord progression.
    """
    song_title = "Michelle"
    
    # The song starts on the tonic chord in the key of F.
    tonic_chord = "F Major"
    
    # It then "jumps" to the tonic minor chord. This I -> i (major to minor tonic)
    # change is very distinctive and likely what is being described.
    next_chord = "F Minor"
    
    print(f"The song by The Beatles that famously starts with a chord 'jump' from the tonic major to the tonic minor is '{song_title}'.")
    print("\nThis chord change (from I to i) is likely what the user means by 'jump from the tonic to the minor fifth chord'.")
    print("\nThe specific chord change at the beginning of the song is:")
    
    # As requested, printing each part of the final progression "equation".
    print(tonic_chord)
    print("->")
    print(next_chord)

find_beatles_song()