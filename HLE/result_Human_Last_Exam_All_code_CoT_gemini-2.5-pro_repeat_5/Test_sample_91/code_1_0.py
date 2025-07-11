def explain_bebop_scale():
    """
    Provides an explanation of the Bebop Scale and prints an example.
    This is the type of scale Erroll Garner uses in his virtuosic runs.
    """
    scale_type = "Bebop Major Scale"
    explanation = (
        "The scale used in the right-hand melody is a Bebop Scale.\n"
        "This type of scale adds a chromatic passing note to a standard 7-note scale, "
        "creating an 8-note scale that aligns well with jazz rhythms.\n"
    )

    # Example: C Major Bebop scale
    # It's a C Major scale (C, D, E, F, G, A, B) with an added G# (a chromatic passing tone).
    root_note = "C"
    notes = ["C", "D", "E", "F", "G", "G#", "A", "B"]

    print(explanation)
    print(f"An example is the {root_note} {scale_type}.")
    print("The 'equation' of notes in this scale is:")
    
    # The prompt asks to "output each number in the final equation".
    # We will represent this by printing the sequence of notes.
    final_equation = " + ".join(notes)
    print(final_equation)

explain_bebop_scale()