def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """

    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # It consists of Jins Ajam on the root, Jins Hijaz on the 4th, and a final whole tone to the octave.
    # Jins Ajam intervals: Whole (1), Whole (1), Half (0.5)
    # Jins Hijaz intervals: Half (0.5), Whole-and-a-half (1.5), Half (0.5)
    # Final connecting interval: Whole (1)
    ascending_intervals = [1, 1, 0.5, 0.5, 1.5, 0.5, 1]

    # Step 2: Define the intervals for the modified descent.
    # The upper register is replaced with Jins Nahawand on the 4th.
    # The notes of this jins, starting on the 4th degree (F), would be F, G, Ab, Bb.
    # The full upper part of the scale becomes F, G, Ab, Bb, C'.
    # The descent is from C' to F: C' -> Bb -> Ab -> G -> F.
    # Interval C' -> Bb: Whole (1)
    # Interval Bb -> Ab: Whole (1)
    # Interval Ab -> G: Half (0.5)
    # Interval G -> F: Whole (1)
    descending_intervals = [1, 1, 0.5, 1]

    # Step 3: Combine the ascending and descending intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output as requested.
    # The list is converted to a string with comma-separated values inside curly braces.
    # The problem states "output each number in the final equation", which means printing the final sequence.
    # All values are already rounded to the nearest quarter interval as required.
    
    # Create the string representation of the list
    interval_str_list = [str(i) for i in all_intervals]
    
    # Format the final output string
    final_output = "{" + ",".join(interval_str_list) + "}"
    
    print(final_output)

solve_maqam_intervals()
<<< {1,1,0.5,0.5,1.5,0.5,1,1,1,0.5,1}>>>