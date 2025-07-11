def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals as described.
    """
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # It consists of Jins Ajam on C (W, W, H), a whole step link, and Jins Hijaz on G (H, A2, H).
    # W = 1, H = 0.5, A2 = 1.5
    # C->D(1), D->E(1), E->F(0.5), F->G(1), G->Ab(0.5), Ab->B(1.5), B->C'(0.5)
    ascending_intervals = [1, 1, 0.5, 1, 0.5, 1.5, 0.5]

    # Step 2: Define the intervals for the descending part.
    # The descent uses a modified scale with Jins Nahawand on F.
    # Jins Nahawand on F gives the notes F, G, Ab, Bb (intervals W, H, W).
    # The descent path is from C' down to F: C' -> Bb -> Ab -> G -> F.
    # C'->Bb (1), Bb->Ab (1), Ab->G (0.5), G->F (1)
    descending_intervals = [1, 1, 0.5, 1]

    # Step 3: Combine the lists to get the full sequence of 11 intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output string as {n1, n2, ...}
    # The str() function is used to convert each number to a string.
    # The ", ".join() method creates a single string with elements separated by ", ".
    result_string = "{" + ", ".join(str(interval) for interval in all_intervals) + "}"
    
    print(result_string)

solve_maqam_intervals()
<<<{"1, 1, 0.5, 1, 0.5, 1.5, 0.5, 1, 1, 0.5, 1"}>>>