def find_beatles_song():
    """
    Identifies and explains the Beatles song with the distinctive chord change at the beginning.
    The user asked for a tonic to a minor fifth (I -> v).
    The actual famous progression is a tonic to a minor fourth (I -> iv).
    """
    song_title = "Michelle"
    key = "F Major"
    tonic_chord = "F Major"
    unusual_chord = "B-flat Minor"
    
    progression_in_numerals = "I -> iv"
    tonic_numeral = 1
    unusual_chord_numeral = 4

    print(f"The song you are likely thinking of is '{song_title}' by The Beatles.")
    print("\nWhile a jump from a tonic to a minor fifth (I -> v) is rare,")
    print("this song features a very famous and similar-sounding jump from a")
    print("tonic major chord to a minor fourth chord.\n")

    print("Here is the analysis:")
    print(f"  - Key: {key}")
    print(f"  - First chord (Tonic, I): {tonic_chord}")
    print(f"  - Second chord (Minor Fourth, iv): {unusual_chord}\n")
    
    # Per the instruction: "output each number in the final equation!"
    # The 'equation' here is the Roman numeral progression.
    print(f"The chord progression can be written with Roman numerals as:")
    print(f"I -> iv")
    print("Breaking down the numerals from that progression:")
    print(f"The first chord is number: {tonic_numeral}")
    print(f"The second chord is number: {unusual_chord_numeral}")


if __name__ == '__main__':
    find_beatles_song()