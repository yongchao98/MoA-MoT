def calculate_maqam_intervals():
    """
    Calculates the sequence of musical intervals based on the musician's performance
    of Maqam Zanjaran, ascending and descending.

    The performance consists of two parts:
    1. Ascending Maqam Zanjaran: 8 notes (C to C'), creating 7 intervals.
       - The scale is Jins Zanjaran on C (C, D-flat, E, F) and Jins Nahawand on F (F, G, A-flat, B-flat).
       - Ascending intervals are calculated based on this structure.
    2. Descending Phrase: 4 notes sung after the ascent, ending on the 4th degree (F).
       - This implies the path from the octave (C') to the 4th degree (F): C' -> B-flat -> A-flat -> G -> F.
       - This adds 4 more intervals to the sequence.

    All interval sizes are recorded as positive values rounded to the nearest 0.25.
    """

    # 1. Calculate the 7 ascending intervals for Maqam Zanjaran (C to C')
    # Jins Zanjaran (Hijaz) on C: C -> Db (0.5), Db -> E (1.5), E -> F (0.5)
    # Jins Nahawand on F: F -> G (1.0), G -> Ab (0.5), Ab -> Bb (1.0)
    # Final step to octave: Bb -> C' (1.0)
    ascending_intervals = [
        0.5,  # C to D-flat (Semitone)
        1.5,  # D-flat to E (Augmented Second, 1.5 tones)
        0.5,  # E to F (Semitone)
        1.0,  # F to G (Whole tone)
        0.5,  # G to A-flat (Semitone)
        1.0,  # A-flat to B-flat (Whole tone)
        1.0   # B-flat to C' (Whole tone)
    ]

    # 2. Calculate the 4 descending intervals
    # The path is from C' down to F: C' -> B-flat -> A-flat -> G -> F
    # The interval sizes are based on the upper part of the scale (Jins Nahawand)
    descending_intervals = [
        1.0,  # C' down to B-flat (Whole tone)
        1.0,  # B-flat down to A-flat (Whole tone)
        0.5,  # A-flat down to G (Semitone)
        1.0   # G down to F (Whole tone)
    ]

    # 3. Combine the intervals into the final sequence of 11 intervals
    total_intervals = ascending_intervals + descending_intervals

    # 4. Format the result as requested
    # The format is '{n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11}'
    # We use a list comprehension and string join to create the final string.
    formatted_output = "{" + ",".join(map(str, total_intervals)) + "}"
    
    print(formatted_output)

calculate_maqam_intervals()
<<<0.5,1.5,0.5,1.0,0.5,1.0,1.0,1.0,1.0,0.5,1.0>>>