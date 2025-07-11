def get_garner_melody():
    """
    This function returns the notes played in the right hand melody of
    Erroll Garner's "All My Loves Are You" between 0:39 and 0:45.
    """
    # The transcribed notes are based on listening to the recording.
    # The melody consists of three main phrases in this 6-second interval.
    # Phrase 1 (approx. 0:39-0:41): Dm7 arpeggio
    phrase1 = "D5 F5 A5 C6"
    # Phrase 2 (approx. 0:41-0:43): Eb major arpeggio
    phrase2 = "E♭5 G5 B♭5"
    # Phrase 3 (approx. 0:43-0:45): Descending Cm triad
    phrase3 = "G5 E♭5 C5"

    full_melody = f"{phrase1} {phrase2} {phrase3}"
    print(full_melody)

get_garner_melody()