import math

def solve():
    """
    Calculates and prints the sequence of musical intervals.
    """
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # The scale is composed of Jins Hijaz on the root (C) and Jins Nahawand on the fifth (G).
    # Jins Hijaz (C-Db-E-F) has intervals: 0.5 (semitone), 1.5 (augmented second), 0.5 (semitone).
    jins_hijaz_intervals = [0.5, 1.5, 0.5]

    # The linking interval between the first jins (ending on F) and the second (starting on G) is a whole tone.
    linking_interval = 1.0

    # Jins Nahawand on G (G-A-Bb-C) has intervals: 1.0 (whole tone), 0.5 (semitone), 1.0 (whole tone).
    jins_nahawand_intervals = [1.0, 0.5, 1.0]

    # The 7 ascending intervals are the combination of these parts.
    ascending_intervals = jins_hijaz_intervals + [linking_interval] + jins_nahawand_intervals
    
    # Step 2: Define the intervals for the descent.
    # The descent uses a modified scale with Jins Nahawand on the fourth degree (F).
    # This implies the upper notes for the descent are F, G, Ab, Bb, C'.
    # The path is from the octave (C') down to the fourth degree (F), singing C' -> Bb -> Ab -> G -> F.
    # The sizes of these intervals (always positive) are:
    # C' to Bb: whole tone (1.0)
    # Bb to Ab: semitone (0.5)
    # Ab to G: semitone (0.5)
    # G to F: whole tone (1.0)
    descending_intervals = [1.0, 0.5, 0.5, 1.0]

    # Step 3: Combine the ascending and descending intervals for the total sequence.
    total_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the list into the required string format '{n1,n2,...}' and print it.
    # The values are already multiples of 0.25, so no rounding is needed.
    # Each number is individually converted to a string and then joined.
    interval_strings = [str(num) for num in total_intervals]
    final_output = f"{{{','.join(interval_strings)}}}"
    
    print(final_output)

solve()
<<<{0.5,1.5,0.5,1.0,1.0,0.5,1.0,1.0,0.5,0.5,1.0}>>>