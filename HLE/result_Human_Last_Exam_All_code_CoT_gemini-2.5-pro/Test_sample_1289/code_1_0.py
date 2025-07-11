def calculate_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # It consists of Jins Ajam on the tonic (1, 1, 0.5), Jins Nahawand on the 4th (1, 0.5, 1),
    # and a final whole tone to the octave (1).
    ascending_intervals = [1, 1, 0.5, 1, 0.5, 1, 1]

    # Step 2: Define the intervals for the descent from the octave to the 4th degree.
    # The notes are C' -> Bb -> Ab -> G -> F.
    # Intervals are (C'-Bb)=1, (Bb-Ab)=1, (Ab-G)=0.5, (G-F)=1.
    descending_intervals = [1, 1, 0.5, 1]

    # Step 3: Combine the ascending and descending intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the list into the required string format.
    # The map(str, ...) converts each number in the list to its string representation.
    # ','.join(...) then concatenates these strings with a comma in between.
    result_string = "{" + ",".join(map(str, all_intervals)) + "}"

    print(result_string)

calculate_maqam_intervals()
print("<<<{" + ",".join(map(str, [1, 1, 0.5, 1, 0.5, 1, 1, 1, 1, 0.5, 1])) + "}>>>")