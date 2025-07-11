def solve_music_intervals():
    """
    This function calculates and prints the sequence of musical intervals
    as described in the problem.
    """

    # 1. Define the 7 intervals of the ascending Maqam Zanjaran.
    # Maqam Zanjaran consists of Jins Rast on the tonic and Jins Hijaz on the 5th.
    # The intervals are [Tone, 3/4 Tone, 3/4 Tone] + [Linking Tone] + [Semitone, 1.5 Tones, Semitone].
    ascending_intervals = [
        1.0,    # Interval 1: 1st to 2nd note (Tone)
        0.75,   # Interval 2: 2nd to 3rd note (3/4 Tone)
        0.75,   # Interval 3: 3rd to 4th note (3/4 Tone)
        1.0,    # Interval 4: 4th to 5th note (Linking Tone)
        0.5,    # Interval 5: 5th to 6th note (Semitone)
        1.5,    # Interval 6: 6th to 7th note (Tone and a half)
        0.5     # Interval 7: 7th to 8th note (Semitone)
    ]

    # 2. Define the 4 intervals of the descending part.
    # The descent uses a scale with Jins Nahawand on the 4th degree.
    # The scale's upper notes are: 4th, 5th, 6th, 7th, 8th(octave).
    # The intervals between these notes are:
    # 4th<->5th (Tone), 5th<->6th (Semitone), 6th<->7th (Tone). This defines the Jins Nahawand.
    # The interval from the 7th to the 8th(octave) in this context is a Tone.
    # The musician descends from the 8th note to the 4th note.
    # Interval sizes are always positive.
    descending_intervals = [
        1.0,    # Interval 8: from 8th to 7th note (Tone)
        1.0,    # Interval 9: from 7th to 6th note (Tone)
        0.5,    # Interval 10: from 6th to 5th note (Semitone)
        1.0     # Interval 11: from 5th to 4th note (Tone)
    ]

    # 3. Combine the lists to get all 11 intervals in order.
    all_intervals = ascending_intervals + descending_intervals

    # 4. Format the output string as specified.
    # The f-string formatting `{x:.2f}` handles cases like 1.0 properly.
    # We then strip any trailing ".0" or useless zeros to get clean numbers.
    formatted_intervals = [f"{x:.2f}".rstrip('0').rstrip('.') for x in all_intervals]
    
    # Create the final string in the format {n1,n2,n3,...}
    final_output_string = "{" + ",".join(formatted_intervals) + "}"
    
    # Print the result to the console.
    print(final_output_string)

solve_music_intervals()
<<<{"1,0.75,0.75,1,0.5,1.5,0.5,1,1,0.5,1"}>>>