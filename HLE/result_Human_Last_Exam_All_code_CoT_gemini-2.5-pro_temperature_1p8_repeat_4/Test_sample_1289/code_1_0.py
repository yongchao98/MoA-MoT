def solve_music_intervals():
    """
    Calculates and prints the sequence of musical intervals sung by the musician.
    """
    # Step 1 & 2: Define and calculate ascending intervals for Maqam Zanjaran.
    # Maqam Zanjaran = Jins Hijaz on root + Jins Nahawand on the 4th.
    # Jins Hijaz (root): 3/4 tone, 3/4 tone, 1/2 tone
    # Connecting interval (4th to 5th degree): whole tone
    # Jins Nahawand (on 4th, starts from 5th): 1/2 tone, whole tone, whole tone
    # This structure is slightly different from some interpretations, but fits the common model.
    # Let's re-evaluate based on the common scale notes: C, D(half-flat), E, F, G, Ab, Bb, C'.
    # 1. C to D(hf): 0.75 (three-quarter tone)
    # 2. D(hf) to E: 0.75 (three-quarter tone)
    # 3. E to F: 0.5 (semitone)
    # 4. F to G: 1.0 (whole tone - link)
    # 5. G to Ab: 0.5 (semitone)
    # 6. Ab to Bb: 1.0 (whole tone)
    # 7. Bb to C': 1.0 (whole tone)
    
    ascending_intervals = [0.75, 0.75, 0.5, 1.0, 0.5, 1.0, 1.0]

    # Step 3: Determine the descending intervals.
    # The musician descends from the octave (note 8) down to the 4th degree (note 4, which is F).
    # The notes of the descent are C', Bb, Ab, G, F. This creates 4 intervals.
    # The scale used is the same, so the intervals are the reverse of the last 4 ascending intervals.
    # Last 4 ascending intervals correspond to the path F -> G -> Ab -> Bb -> C'.
    # These intervals are: [1.0, 0.5, 1.0, 1.0]
    # The descending intervals (C'->Bb, Bb->Ab, Ab->G, G->F) have the same sizes, just in reverse sequence of notes.
    # Interval C' -> Bb has the same size as Bb -> C' (1.0)
    # Interval Bb -> Ab has the same size as Ab -> Bb (1.0)
    # Interval Ab -> G has the same size as G -> Ab (0.5)
    # Interval G -> F has the same size as F -> G (1.0)
    
    descending_intervals = [1.0, 1.0, 0.5, 1.0]
    
    # Step 4: Combine the lists and print in the specified format.
    all_intervals = ascending_intervals + descending_intervals
    
    # Format the list as a string like '{1,2,3}'
    # The problem statement requires formatting with the numbers directly.
    # So we format the final list of numbers.
    # The list is: 0.75, 0.75, 0.5, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.5, 1.0
    
    # Convert numbers to strings to handle .5 and .75 without extra zeros
    formatted_intervals = []
    for num in all_intervals:
        if num == int(num):
            formatted_intervals.append(str(int(num)))
        else:
            formatted_intervals.append(str(num))
            
    result_string = "{" + ",".join(formatted_intervals) + "}"
    
    print("The sequence of intervals is:")
    print(result_string)

solve_music_intervals()

# The final answer in the required format
final_answer_list = [0.75, 0.75, 0.5, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.5, 1.0]
formatted_final_answer = "{" + ",".join(map(str, final_answer_list)) + "}"
print(f"\n<<<{formatted_final_answer}>>>")