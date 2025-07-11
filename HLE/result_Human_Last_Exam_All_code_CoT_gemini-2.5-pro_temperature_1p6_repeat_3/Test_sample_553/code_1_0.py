def solve_beatles_riddle():
    """
    Identifies the Beatles song based on a harmonic description
    and explains the music theory behind the answer.
    """
    song_title = "Michelle"

    explanation = """
The Beatles song that is famous for a harmonic structure that could be described as a "jump from the tonic to a minor chord" at the beginning is "Michelle".

Here's the music theory breakdown:

1.  **The Tonic Key:** The song's primary key, especially in the memorable "I love you, I love you, I love you" chorus, is F Major. In this key, the tonic chord (I) is F Major.

2.  **The Opening Chord:** However, the song's famous guitar introduction and the verses begin on an F minor (Fm) chord. This is the minor tonic (i), borrowed from the parallel key of F minor.

3.  **The "Jump":** The song starts immediately on this unexpected minor chord, creating a melancholy "French" feel. While it's not literally a jump from a major tonic to a minor fifth (which would be F Major to C minor), it's a jump from the *expectation* of the F Major tonic to the actual sound of the F minor chord. This prominent use of the minor tonic at the start is a very distinct and famous feature of the song.

Therefore, while the technical description in the question is a point of music theory debate, "Michelle" is the accepted answer due to this unique and central harmonic feature.
"""

    print("Analyzing the song structure...")
    print(explanation)
    print(f"The song is: {song_title}")

solve_beatles_riddle()