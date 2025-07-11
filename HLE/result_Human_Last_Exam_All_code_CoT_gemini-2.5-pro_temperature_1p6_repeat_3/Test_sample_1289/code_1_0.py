def solve_music_intervals():
    """
    This function calculates and prints the sequence of musical intervals
    as described by the user's prompt.
    """

    # Define the interval values for clarity
    semitone = 0.5
    three_quarter_tone = 0.75
    whole_tone = 1.0

    # Step 1: Define the 7 ascending intervals for Maqam Zanjaran.
    # The structure is Jins Rast, Whole Tone, Jins Rast.
    # Jins Rast = [Whole, Three-Quarter, Three-Quarter]
    ascending_intervals = [
        whole_tone,          # Interval 1
        three_quarter_tone,  # Interval 2
        three_quarter_tone,  # Interval 3
        whole_tone,          # Interval 4 (connecting)
        whole_tone,          # Interval 5 (from Jins Rast on the 5th)
        three_quarter_tone,  # Interval 6
        three_quarter_tone,  # Interval 7
    ]

    # Step 2: Define the 4 descending intervals using the modified scale.
    # The descent is from the 8th to the 4th note, using a scale with
    # Jins Nahawand on the 4th. The path is C' -> Bb -> Ab -> G -> F.
    descending_intervals = [
        whole_tone,          # Interval 8: C' to Bb
        whole_tone,          # Interval 9: Bb to Ab
        semitone,            # Interval 10: Ab to G
        whole_tone,          # Interval 11: G to F
    ]

    # Step 3: Combine the lists to get all 11 intervals in order.
    all_intervals = ascending_intervals + descending_intervals

    # Helper function to format numbers cleanly (e.g., 1.0 becomes "1")
    def format_number(n):
        if n == int(n):
            return str(int(n))
        return str(n)

    # Format the final list into the specified string format: "{n1,n2,n3,...}"
    # This also fulfills the requirement to "output each number in the final equation".
    output_string = "{" + ",".join(map(format_number, all_intervals)) + "}"

    print(output_string)

solve_music_intervals()
<<<{"1,0.75,0.75,1,1,0.75,0.75,1,1,0.5,1"}>>>