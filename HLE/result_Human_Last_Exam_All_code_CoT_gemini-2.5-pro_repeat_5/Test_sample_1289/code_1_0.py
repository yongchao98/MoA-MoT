def solve_music_intervals():
    """
    Calculates and prints the sequence of musical intervals as described.
    """

    # Step 1 & 2: Define the intervals based on the musical scales.

    # The ascending phrase uses Maqam Zanjaran.
    # It consists of Jins Rast on the tonic (C -> F) and Jins Hijaz on the 5th (G -> C).
    # Jins Rast (C,D,E-half-flat,F) intervals: Whole (1), Three-quarter (0.75), Three-quarter (0.75)
    # Connecting interval (F -> G): Whole (1)
    # Jins Hijaz (G,A-flat,B,C) intervals: Half (0.5), One-and-a-half (1.5), Half (0.5)
    ascending_intervals = [1, 0.75, 0.75, 1, 0.5, 1.5, 0.5]

    # The descending phrase uses a scale with Jins Nahawand on the 4th degree (F).
    # The musician sings 8 notes up to the octave (C), then sings 4 notes down,
    # ending on the 4th degree (F). The four notes of the descent are B-flat, A-flat, G, F.
    
    # The intervals for the descent are calculated as follows:
    # 1. Transition from octave C down to B-flat: Whole tone (1)
    # 2. B-flat down to A-flat: Whole tone (1)
    # 3. A-flat down to G: Half tone (0.5)
    # 4. G down to F: Whole tone (1)
    descending_intervals = [1, 1, 0.5, 1]

    # Step 3: Combine all intervals into a single list.
    # There are 7 ascending intervals and 4 intervals in the descending part (1 transition + 3 descent).
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format and print the result as specified.
    # The problem asks to output each number in the final list.
    interval_strings = []
    for interval in all_intervals:
      # Format to remove trailing .0 for whole numbers like '1.0'
      if interval == int(interval):
        interval_strings.append(str(int(interval)))
      else:
        interval_strings.append(str(interval))
        
    final_output = "{" + ",".join(interval_strings) + "}"
    
    print(final_output)

solve_music_intervals()