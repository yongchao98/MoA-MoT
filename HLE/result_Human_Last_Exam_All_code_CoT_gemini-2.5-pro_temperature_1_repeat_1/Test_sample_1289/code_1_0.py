def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals
    as described in the problem.
    """

    # Step 1: Define the 7 ascending intervals for Maqam Zanjaran (Major Scale).
    # Pattern: Whole, Whole, Half, Whole, Whole, Whole, Half
    ascending_intervals = [1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.5]

    # Step 2: Define the 4 descending intervals using the modified scale.
    # The descent is from the octave (8) to the 4th degree, using a modified
    # upper register based on Jins Nahawand on the 4th degree.
    # This results in the following interval pattern for the descent from 8->7'->6'->5->4:
    # 8 -> 7' (e.g., C -> Bb) = Whole Tone = 1.0
    # 7' -> 6' (e.g., Bb -> Ab) = Whole Tone = 1.0
    # 6' -> 5  (e.g., Ab -> G) = Semi-tone = 0.5
    # 5 -> 4   (e.g., G -> F) = Whole Tone = 1.0
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the ascending and descending intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the list for the final output string.
    # Convert floats to integers where possible for cleaner output (e.g., 1.0 -> 1).
    formatted_parts = []
    for interval in all_intervals:
        if interval == int(interval):
            formatted_parts.append(str(int(interval)))
        else:
            formatted_parts.append(str(interval))

    # The final format is {n1,n2,n3,...}
    output_string = "{" + ",".join(formatted_parts) + "}"

    # Print the final result.
    print(output_string)

solve_maqam_intervals()
<<< {1,1,0.5,1,1,1,0.5,1,1,0.5,1} >>>