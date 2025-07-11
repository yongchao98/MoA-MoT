def solve_music_theory_task():
    """
    This script identifies and displays the musical scale used in a specific
    part of an Erroll Garner performance.
    """
    
    # The scale is identified by analyzing the notes in the phrase.
    # The tonal center is C.
    # The fast descending run from 0:39-0:40 uses the notes of the C Bebop Dominant scale.
    # This scale is a C Mixolydian scale (C D E F G A Bb) with an added chromatic
    # passing tone, which is the major 7th (B).
    scale_name = "C Bebop Dominant Scale"
    
    # The notes of the C Bebop Dominant Scale
    # In music, these notes correspond to the scale degrees: 1, 2, 3, 4, 5, 6, b7, 7
    notes = ["C", "D", "E", "F", "G", "A", "Bb", "B"]
    
    print(f"The scale Erroll Garner plays in the right hand melody is the {scale_name}.")
    print("This scale adds a chromatic passing tone to the standard Mixolydian mode,")
    print("which is a hallmark of the bebop jazz style.")
    print("\nThe notes in the scale are:")
    
    # Print each note of the scale as requested
    output_string = " ".join(notes)
    print(output_string)

solve_music_theory_task()