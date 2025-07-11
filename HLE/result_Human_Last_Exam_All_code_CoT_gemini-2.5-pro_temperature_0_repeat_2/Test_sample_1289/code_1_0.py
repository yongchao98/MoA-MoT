def solve_maqam_intervals():
    """
    This function calculates and prints the sequence of musical intervals
    based on the provided description.
    """
    # Step 1: Define the 7 ascending intervals of Maqam Zanjaran.
    # The scale structure is Jins Ajam on C + Jins Hijaz on G.
    # Notes: C, D, E, F, G, Ab, B, C'
    # Intervals:
    # C-D: Whole tone (1.0)
    # D-E: Whole tone (1.0)
    # E-F: Semitone (0.5)
    # F-G: Whole tone (1.0)
    # G-Ab: Semitone (0.5)
    # Ab-B: Augmented second (1.5)
    # B-C': Semitone (0.5)
    ascending_intervals = [1.0, 1.0, 0.5, 1.0, 0.5, 1.5, 0.5]

    # Step 2: Define the 4 descending intervals.
    # The descent uses a scale with Jins Nahawand on the 4th degree (F).
    # The upper notes of this scale are F, G, Ab, Bb, C'.
    # The musician descends from the 8th note (C') to the 4th (F).
    # Notes sung: C', Bb, Ab, G, F
    # Intervals:
    # C' down to Bb: Whole tone (1.0)
    # Bb down to Ab: Whole tone (1.0)
    # Ab down to G: Semitone (0.5)
    # G down to F: Whole tone (1.0)
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the intervals into one list.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the list into the required string format.
    # We convert numbers to integers if they have no fractional part for cleaner output.
    formatted_numbers = []
    for num in all_intervals:
        if num == int(num):
            formatted_numbers.append(str(int(num)))
        else:
            formatted_numbers.append(str(num))

    # Join the numbers with commas and enclose in curly braces.
    output_string = "{" + ",".join(formatted_numbers) + "}"

    print(output_string)

solve_maqam_intervals()
<<<{"1,1,0.5,1,0.5,1.5,0.5,1,1,0.5,1"}>>>