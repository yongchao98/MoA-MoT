def solve_intervals():
    """
    Calculates and prints the sequence of musical intervals as described.
    """

    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # Jins `Ajam (1, 1, 0.5) + connecting tone (1) + Jins Nahawand on the 5th (1, 0.5, 1)
    ascending_intervals = [1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0]

    # Step 2: Define the intervals for the descending modified scale.
    # The descent is from the octave to the 4th degree, using a scale
    # with Jins Nahawand on the 4th degree.
    # Path: Octave(C') -> Bb -> Ab -> G -> F(4th)
    # Intervals: C'-Bb (1), Bb-Ab (1), Ab-G (0.5), G-F (1)
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the two lists to get the full sequence of 11 intervals.
    total_intervals = ascending_intervals + descending_intervals

    # Helper function to format numbers as integers if they are whole, else as floats.
    def format_number(n):
        if n == int(n):
            return str(int(n))
        return str(n)

    # Format the final list into the required string format "{n1,n2,...}".
    formatted_intervals = ",".join(map(format_number, total_intervals))
    final_output = f"{{{formatted_intervals}}}"

    print(final_output)

solve_intervals()
