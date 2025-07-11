def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for the specified piano chord.
    """
    # The chord notes are Bb3, Db4, F4, Bb4, Db5.
    # The span is a minor tenth, which is very wide.
    # 1 (thumb) on the lowest note (Bb3).
    # 5 (pinky) on the highest note (Db5).
    # 2, 3, 4 are assigned sequentially to the inner notes.
    fingering = "Bb3 1, Db4 2, F4 3, Bb4 4, Db5 5"
    print(fingering)

get_piano_fingering()