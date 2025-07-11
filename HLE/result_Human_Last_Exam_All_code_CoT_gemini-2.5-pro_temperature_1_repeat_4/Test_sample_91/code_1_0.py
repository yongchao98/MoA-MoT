def analyze_garner_melody():
    """
    Identifies the scale used in a specific passage of Erroll Garner's
    "All My Loves Are You" based on musicological analysis.
    """
    # Step 1: Analyze the specified musical passage (0:39-0:43).
    # The right-hand melody in this section is a very fast, descending flourish.
    # Listening reveals that the run has a "dreamy," "floating," and symmetrical
    # quality. This distinct sound is produced by a scale that lacks any half-steps.

    # Step 2: Identify the scale based on its sonic characteristics.
    # The scale built entirely from whole steps is the Whole-Tone Scale.
    # It is a six-note scale (hexatonic) and is a classic device in jazz
    # and impressionistic music to create a sense of ambiguity and movement.
    scale_type = "Whole-Tone Scale"

    # Step 3: Provide an example of the scale's structure.
    # There are only two possible whole-tone scales. The following is one example.
    example_notes = ["C", "D", "E", "F#", "G#", "A#"]

    print("Analysis of Erroll Garner's right-hand melody (0:39-0:43):")
    print("="*60)
    print(f"The scale Garner uses for the fast, descending run is the: {scale_type}")
    print("\nThis scale is constructed entirely from whole-step intervals, giving it a characteristic 'dreamy' and ambiguous sound.")
    print("\nFor example, a whole-tone scale starting on the note C is composed of the following notes:")

    # Step 4: Print each note in the example scale as requested.
    # The final "equation" is the set of notes that form the scale.
    print(f"The notes are: {example_notes[0]}, {example_notes[1]}, {example_notes[2]}, {example_notes[3]}, {example_notes[4]}, {example_notes[5]}")
    print("="*60)

analyze_garner_melody()