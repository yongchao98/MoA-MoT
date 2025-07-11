def find_lowest_note():
    """
    This function identifies the lowest of the final three notes played by the
    first violin in Dvořák's Symphony No. 2.
    """
    # According to the musical score, the first violin part ends on a
    # three-note chord (a triple-stop).
    note1 = "D"
    note2 = "F"
    note3 = "B-flat"

    # The notes are played simultaneously. We need to identify the lowest in pitch.
    lowest_note = "D"

    print(f"The final three notes played by the first violin are {note1}, {note2}, and {note3}.")
    print(f"The lowest of these notes is: {lowest_note}")

find_lowest_note()