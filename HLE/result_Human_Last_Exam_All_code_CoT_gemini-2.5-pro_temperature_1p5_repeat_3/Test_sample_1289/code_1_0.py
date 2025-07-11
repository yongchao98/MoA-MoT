def solve_maqam_intervals():
    """
    This function calculates and prints the sequence of musical intervals
    as described in the problem.
    """
    
    # Step 1: Determine the 7 ascending intervals of Maqam Zanjaran.
    # The scale is composed of Jins Rast on the root (C) and Jins Hijaz on the 5th (G).
    # Jins Rast on C gives notes C, D, E-half-flat, F.
    # Intervals: C-D (Whole=1), D-E-hf (3/4=0.75), E-hf-F (3/4=0.75).
    # The linking interval F-G is a Whole tone (1).
    # Jins Hijaz on G gives notes G, A-flat, B, C'.
    # Intervals: G-Ab (Semi=0.5), Ab-B (1.5), B-C' (Semi=0.5).
    ascending_intervals = [1, 0.75, 0.75, 1, 0.5, 1.5, 0.5]

    # Step 2: Determine the 4 descending intervals.
    # The descent is from the octave (C') to the fourth degree (F).
    # The scale modification is a Jins Nahawand on the 4th degree (F).
    # Jins Nahawand on F gives notes F, G, Ab, Bb. Intervals: F-G(1), G-Ab(0.5), Ab-Bb(1).
    # This implies the upper notes of the descending scale are based on this.
    # To complete the scale to the octave C', the interval from Bb is a whole tone (1).
    # The notes sung in descent are: C', Bb, Ab, G, F.
    # The intervals are calculated for this descending sequence:
    # Interval C' -> Bb is a whole tone (1).
    # Interval Bb -> Ab is a whole tone (1).
    # Interval Ab -> G is a semitone (0.5).
    # Interval G -> F is a whole tone (1), from the unmodified lower part of the scale.
    descending_intervals = [1, 1, 0.5, 1]
    
    # Step 3: Combine the lists and format the output.
    all_intervals = ascending_intervals + descending_intervals

    # Format the list into the required string format: "{n1,n2,...}"
    # The problem statement requires printing each number in the final equation.
    # We will format the list of floats into a comma-separated string.
    # The 'g' format specifier is used to avoid trailing zeros for numbers like 1.0.
    formatted_intervals = ",".join(f"{n:g}" for n in all_intervals)
    
    output_string = f"{{{formatted_intervals}}}"
    
    print(output_string)

solve_maqam_intervals()
<<<1,0.75,0.75,1,0.5,1.5,0.5,1,1,0.5,1>>>