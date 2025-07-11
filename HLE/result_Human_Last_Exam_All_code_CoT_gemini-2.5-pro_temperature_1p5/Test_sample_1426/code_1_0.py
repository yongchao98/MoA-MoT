def find_lowest_note():
    """
    Identifies the final three notes for the first violin in Dvořák's 
    Symphony No. 2, Op. 4, and determines the lowest one.
    """
    
    # The final three notes played as a chord by the first violins.
    # In scientific pitch notation, they are B-flat 4, D 5, and B-flat 5.
    final_notes = ["B-flat", "D", "B-flat (octave higher)"]
    lowest_note = "B-flat"

    print("At the end of Dvořák's Symphony No. 2, Op. 4, the first violin section plays a final chord.")
    print("The notes in this final chord are:")
    for note in final_notes:
        print(f"- {note}")
    
    print("\nComparing these pitches, the lowest note is:")
    print(lowest_note)

find_lowest_note()