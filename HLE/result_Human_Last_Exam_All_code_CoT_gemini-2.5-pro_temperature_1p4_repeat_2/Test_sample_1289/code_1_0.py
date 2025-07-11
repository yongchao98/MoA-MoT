def solve_musical_intervals():
    """
    Calculates and prints a sequence of musical intervals based on Maqam Zanjaran
    and a modified descending scale.
    """
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # It consists of Jins Ajam on the root (1, 1, 0.5) and Jins Hijaz on the 5th (0.5, 1.5, 0.5),
    # linked by a whole tone (1).
    # The seven ascending intervals are:
    # 1. Whole Tone (1.0)
    # 2. Whole Tone (1.0)
    # 3. Semitone (0.5)
    # 4. Whole Tone (1.0) - linking interval
    # 5. Semitone (0.5)
    # 6. Augmented Second (1.5)
    # 7. Semitone (0.5)
    ascending_intervals = [1, 1, 0.5, 1, 0.5, 1.5, 0.5]

    # Step 2: Define the intervals for the descending scale.
    # The upper register is replaced by Jins Nahawand on the 4th degree.
    # Jins Nahawand from the 4th has notes (relative to tonic): F, G, Ab, Bb.
    # The scale reaches the octave C' from Bb (a whole tone).
    # The musician descends from the 8th note (C') to the 4th note (F).
    # The four descending intervals are:
    # 8. C' down to Bb: Whole Tone (1.0)
    # 9. Bb down to Ab: Whole Tone (1.0)
    # 10. Ab down to G: Semitone (0.5)
    # 11. G down to F: Whole Tone (1.0)
    descending_intervals = [1, 1, 0.5, 1]

    # Step 3: Combine the ascending and descending intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Helper function to format numbers (e.g., 1.0 becomes '1', 1.5 remains '1.5')
    def format_number(n):
      if n == int(n):
        return str(int(n))
      return str(n)

    # Step 4: Format the final list into the required string format.
    # Note: The problem asks us to output each number in the final equation.
    # This loop builds the string representation piece by piece.
    result_string = "{"
    for i, interval in enumerate(all_intervals):
        result_string += format_number(interval)
        if i < len(all_intervals) - 1:
            result_string += ", "
    result_string += "}"
    
    print(result_string)

solve_musical_intervals()