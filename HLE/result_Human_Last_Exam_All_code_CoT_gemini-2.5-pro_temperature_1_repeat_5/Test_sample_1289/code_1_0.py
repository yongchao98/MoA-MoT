def solve_music_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the described performance.
    """
    # Step 1: Define the ascending intervals for Maqam Zanjaran.
    # Maqam Zanjaran on C: C D E F G Ab B C
    # Jins Ajam on C (C-D-E-F): Tone, Tone, Semitone -> {1, 1, 0.5}
    # Link (F-G): Tone -> {1}
    # Jins Hijaz on G (G-Ab-B-C): Semitone, Augmented Second, Semitone -> {0.5, 1.5, 0.5}
    ascending_intervals = [1, 1, 0.5, 1, 0.5, 1.5, 0.5]

    # Step 2: Define the descending intervals.
    # The descent is from the octave (C) to the 4th degree (F).
    # The scale uses Jins Nahawand on the 4th degree (F).
    # Jins Nahawand on F (F-G-Ab-Bb): Tone, Semitone, Tone.
    # This defines the upper scale notes for the descent as F, G, Ab, Bb, C.
    # The musician sings descending notes: C, Bb, Ab, G, F.
    # Interval C down to Bb: Tone -> {1}
    # Interval Bb down to Ab: Semitone -> {0.5}
    # Interval Ab down to G: Semitone -> {0.5}
    # Interval G down to F: Tone -> {1}
    descending_intervals = [1, 0.5, 0.5, 1]

    # Step 3: Combine the lists.
    # The total performance has 7 ascending + 4 descending = 11 intervals.
    all_intervals = ascending_intervals + descending_intervals

    # Format the output string as requested, e.g., "{1,2,3}"
    # The map(str, ...) converts each number to a string.
    # ','.join(...) creates the comma-separated string of numbers.
    result_string = "{" + ",".join(map(str, all_intervals)) + "}"
    
    print(result_string)

solve_music_intervals()

# The final answer in the required format
final_answer = "{1,1,0.5,1,0.5,1.5,0.5,1,0.5,0.5,1}"
print(f"<<<{final_answer}>>>")