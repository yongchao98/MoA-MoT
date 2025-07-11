def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals sung by the musician.
    """
    # Step 1 & 2: Define the intervals based on Maqam Zanjaran structure.
    # Jins Hijaz on the root: 0.5, 1.5, 0.5
    # Linking interval: 1.0
    # Jins Nahawand on the 4th: (1.0), 0.5, 1.0
    # The final interval is from the 7th degree to the octave.
    
    # Step 3: Calculate the 7 ascending intervals.
    # Notes: C -> Db -> E -> F -> G -> Ab -> Bb -> C'
    ascending_intervals = [
        0.5,  # C to Db (Semitone)
        1.5,  # Db to E (Augmented Second)
        0.5,  # E to F (Semitone)
        1.0,  # F to G (Whole tone)
        0.5,  # G to Ab (Semitone)
        1.0,  # Ab to Bb (Whole tone)
        1.0   # Bb to C' (Whole tone)
    ]
    
    # Step 4: Calculate the 4 descending intervals.
    # The musician descends from the 8th note to the 4th note (C' -> Bb -> Ab -> G -> F)
    # The intervals are calculated based on the same upper scale structure (Jins Nahawand).
    descending_intervals = [
        1.0,  # C' to Bb (Whole tone)
        1.0,  # Bb to Ab (Whole tone)
        0.5,  # Ab to G (Semitone)
        1.0   # G to F (Whole tone)
    ]
    
    # Step 5: Combine the interval lists.
    all_intervals = ascending_intervals + descending_intervals
    
    # Format the list into a string "{v1,v2,v3,...}"
    # The numbers are already in quarter intervals, so no rounding is needed.
    result_string = "{" + ",".join(map(str, all_intervals)) + "}"
    
    print(result_string)

solve_maqam_intervals()