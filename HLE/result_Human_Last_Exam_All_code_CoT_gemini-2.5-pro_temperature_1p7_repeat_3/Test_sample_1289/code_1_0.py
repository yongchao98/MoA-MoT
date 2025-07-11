def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the user's description.
    """
    # Step 1: Determine the 7 ascending intervals for Maqam Zanjaran.
    # Maqam Zanjaran = Jins Hijaz on tonic + Jins Rast on 5th.
    # Jins Hijaz (e.g., C-Db-E-F) has intervals: semi-tone, augmented second, semi-tone.
    jins_hijaz_intervals = [0.5, 1.5, 0.5]

    # The interval connecting the two jins (e.g., F to G).
    connecting_interval = [1.0]

    # Jins Rast on the 5th (e.g., starting on G). To make the full scale sum to an octave (6 tones),
    # the intervals are: three-quarter tone, whole tone, three-quarter tone.
    jins_rast_intervals = [0.75, 1.0, 0.75]

    ascending_intervals = jins_hijaz_intervals + connecting_interval + jins_rast_intervals

    # Step 2: Determine the 4 descending intervals.
    # The descent is from the octave (C') to the 4th degree (F) using a scale with Jins Nahawand on the 4th degree.
    # This implies a natural minor-like descent in the upper register.
    # The notes are C' -> Bb -> Ab -> G -> F.
    # C' to Bb = whole tone
    # Bb to Ab = whole tone
    # Ab to G = semi-tone
    # G to F = whole tone
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine all 11 intervals in order.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output string as specified.
    # The f-string and join method create the comma-separated list inside the curly braces.
    # The values are already exact quarter-intervals, so no rounding is needed.
    output_string = "{" + ",".join(map(str, all_intervals)) + "}"
    
    print(output_string)

solve_maqam_intervals()