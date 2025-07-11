def calculate_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals as described.
    """
    # Step 1: Define the 7 intervals for the ascending Maqam Zanjaran.
    # This Maqam is composed of Jins Hijaz (0.5, 1.5, 0.5) on the root,
    # a linking interval (1.0), and Jins Nahawand (1.0, 0.5, 1.0) starting on the 5th degree.
    # The full sequence of intervals between the 8 notes is:
    # 1->2: 0.5
    # 2->3: 1.5
    # 3->4: 0.5
    # 4->5: 1.0
    # 5->6: 0.5 (The prompt's description implies Jins Nahawand on the 4th, which means the intervals from the 4th degree are 1, 0.5, 1. So 4->5 is 1, 5->6 is 0.5, 6->7 is 1)
    # Let's correct the sequence based on Jins Nahawand on the 4th.
    # Jins Hijaz on root (notes 1-4): C, Db, E, F -> Intervals: 0.5, 1.5, 0.5
    # Jins Nahawand on 4th (notes 4-8): F, G, Ab, Bb, C -> Intervals: 1, 0.5, 1, 1
    # So the full ascending interval sequence is:
    ascending_intervals = [0.5, 1.5, 0.5, 1.0, 0.5, 1.0, 1.0]

    # Step 2: Define the 4 intervals for the descent from the 8th to the 4th note.
    # The intervals are the reverse of the upper part of the scale.
    # Interval 8->7 (same size as 7->8): 1.0
    # Interval 7->6 (same size as 6->7): 1.0
    # Interval 6->5 (same size as 5->6): 0.5
    # Interval 5->4 (same size as 4->5): 1.0
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the lists to get all 11 intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the list into the required string format.
    # Numbers should be rounded to the nearest quarter, but are already precise.
    # Whole numbers should be printed as integers (e.g., '1' not '1.0').
    formatted_numbers = []
    for interval in all_intervals:
        if interval == int(interval):
            formatted_numbers.append(str(int(interval)))
        else:
            formatted_numbers.append(str(interval))

    final_output = "{" + ",".join(formatted_numbers) + "}"
    print(final_output)

calculate_maqam_intervals()
<<< {0.5,1.5,0.5,1,0.5,1,1,1,1,0.5,1} >>>