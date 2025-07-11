def solve_musical_intervals():
    """
    This function calculates and prints the sequence of musical intervals
    based on the described performance.
    """

    # Step 1: Define the 7 ascending intervals from Maqam Zanjaran.
    # The Maqam is Jins Rast on C (1.0, 0.75, 0.75) linked by a whole tone (1.0)
    # to Jins Hijaz on G (0.5, 1.5, 0.5).
    ascending_intervals = [1.0, 0.75, 0.75, 1.0, 0.5, 1.5, 0.5]

    # Step 2: Define the 4 descending intervals.
    # The descent path is C -> Bb -> Ab -> G -> F.
    # The scale used provides these intervals:
    # C -> Bb (Whole tone)
    # Bb -> Ab (Whole tone)
    # Ab -> G (Semitone)
    # G -> F (Whole tone)
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the two lists into the full sequence of 11 intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output string as '{n1,n2,n3,...}'.
    # We must ensure numbers like 1.0 are printed correctly without trailing zeros
    # where possible, e.g., '1' instead of '1.0'. The map to str handles this.
    # However, the prompt's examples use decimals, so we'll format them consistently.
    interval_strings = [f"{i:.2f}".replace('.00', '.0').replace('.50', '.5').replace('.25', '.25').replace('.75', '.75') for i in all_intervals]
    
    # Final check on formatting based on prompt examples.
    # 0.25, 0.5, 0.75, 1, 1.25 -> My formatting produces '0.25', '0.5', '0.75', '1.0', '1.25'. The prompt example is '1', not '1.0'. Let's adjust.
    final_interval_strings = []
    for i in all_intervals:
        if i == int(i):
            final_interval_strings.append(str(int(i)))
        else:
            final_interval_strings.append(str(i))


    output_string = "{" + ",".join(final_interval_strings) + "}"
    print(output_string)

solve_musical_intervals()
