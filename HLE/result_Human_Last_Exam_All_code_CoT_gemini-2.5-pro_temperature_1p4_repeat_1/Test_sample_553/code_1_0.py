def explain_beatles_harmony():
    """
    Explains which Beatles song starts with a jump from the tonic to the minor fifth chord.
    """
    song_title = "Michelle"
    key = "F minor"
    tonic_chord = "F minor (i)"
    minor_fifth_chord = "C minor (v)"

    explanation = f"""The song by The Beatles that features a distinctive move from the tonic to the minor fifth chord is "{song_title}".

Here is the harmonic breakdown:

1.  The song's verse is in the key of {key}.
2.  The tonic chord (i) is {tonic_chord}.
3.  The minor fifth chord (v) in this key is {minor_fifth_chord}.

The song famously uses the progression from {tonic_chord} to {minor_fifth_chord}, creating its unique, sophisticated sound.
"""
    print(explanation)

explain_beatles_harmony()