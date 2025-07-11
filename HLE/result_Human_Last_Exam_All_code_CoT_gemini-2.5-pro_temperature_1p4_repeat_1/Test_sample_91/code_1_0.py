def solve_music_theory_task():
    """
    Analyzes and identifies the musical scale used in a specific musical passage.

    The task is to identify the scale Erroll Garner plays in the right-hand melody
    of "All My Loves Are You" between seconds 39 and 43.

    Analysis:
    1. Transcription of the melody reveals a descending run.
    2. The notes in this run are: C, B-flat, A, G, F, E-flat, and D.
    3. Organizing these notes starting from the root note C gives us the collection:
       C, D, E-flat, F, G, A, B-flat.
    4. We analyze the intervals from the root C:
       - C to D: Major Second (whole step)
       - C to E-flat: Minor Third (one and a half steps)
       - C to F: Perfect Fourth (two and a half steps)
       - C to G: Perfect Fifth (three and a half steps)
       - C to A: Major Sixth (four and a half steps)
       - C to B-flat: Minor Seventh (five whole steps)
    5. This specific interval pattern (Root, Major 2nd, Minor 3rd, Perfect 4th,
       Perfect 5th, Major 6th, Minor 7th) is the defining characteristic of the
       Dorian mode.
    6. Since the root note is C, the scale is C Dorian.
    """

    scale_name = "C Dorian"
    root_note = "C"
    notes = ["C", "D", "E-flat", "F", "G", "A", "B-flat"]
    # The "equation" can be represented by the scale degrees or interval formula.
    # 1=Root, 2=Major 2nd, b3=Minor 3rd, 4=Perfect 4th, 5=Perfect 5th, 6=Major 6th, b7=Minor 7th
    interval_equation = ["1", "2", "b3", "4", "5", "6", "b7"]
    
    print(f"The scale identified in the passage is: {scale_name}")
    print("-" * 35)
    print("This scale is comprised of the following notes:")
    print("Notes: " + ", ".join(notes))
    print("\nThe structure, or 'equation,' of this scale is defined by its intervals from the root:")
    
    # We print each "number" (interval) in the final equation.
    for i in range(len(notes)):
        print(f"Interval {interval_equation[i]}: {notes[i]}")

solve_music_theory_task()