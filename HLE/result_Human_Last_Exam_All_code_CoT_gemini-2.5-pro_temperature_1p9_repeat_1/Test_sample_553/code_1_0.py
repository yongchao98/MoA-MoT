import textwrap

def solve_beatles_song_query():
    """
    This function identifies and explains which Beatles song fits the
    description of starting with a chord jump from the tonic to the minor fifth.
    """

    # Explanation of the music theory concepts
    explanation = """
    The user is asking for a Beatles song with a specific opening chord progression: from the tonic chord to the minor fifth chord.

    1.  Tonic Chord: This is the 'home' chord of the song's key. For a song in a minor key, this is the 'i' chord.
    2.  Minor Fifth Chord: This is a minor chord built on the fifth note of the key's scale. In a minor key, this is the 'v' chord.
    3.  The Progression: The request describes the chord change i -> v.

    The song that most clearly features this progression near its start is "Michelle".
    - The verses of "Michelle" are in the key of F minor.
    - The tonic chord (i) is F minor (Fm).
    - The minor fifth chord (v) is C minor (Cm).

    In the first verse, during the French lyric "Sont les mots qui vont très bien ensemble", the chord progression is a direct jump from Fm to Cm.
    """

    song_title = "Michelle"

    print("Analyzing the musical query...")
    print(textwrap.dedent(explanation))
    print("This jump from F minor to C minor perfectly matches the 'tonic to the minor fifth chord' description.")
    print("-" * 30)
    print(f"The song is: {song_title}")
    print("-" * 30)


solve_beatles_song_query()

# Final Answer
print("<<<Michelle>>>")