import sys
# A bug in various Python versions will not print large strings to stdout without this
if sys.stdout.encoding != 'UTF-8':
    sys.stdout.reconfigure(encoding='UTF-8')
    
def calculate_maqam_intervals():
    """
    This function calculates the sequence of intervals sung by the musician
    and prints them in the specified format.
    """

    # Step 1: Define the intervals for the ascending scale (Maqam Zanjaran).
    # The scale consists of Jins Hijaz on the root and Jins Nahawand on the 4th.
    
    # Jins Hijaz intervals: semitone, augmented second, semitone
    jins_hijaz = [0.5, 1.5, 0.5]
    
    # Jins Nahawand intervals: whole tone, semitone, whole tone
    jins_nahawand = [1, 0.5, 1]
    
    # The final interval to complete the octave is a whole tone.
    final_ascending_interval = [1.0]
    
    # Combine the parts to get the 7 ascending intervals.
    ascending_intervals = jins_hijaz + jins_nahawand + final_ascending_interval
    
    # Step 2: Define the intervals for the descending scale.
    # The musician descends 4 intervals from the octave to the 4th degree.
    # This path (octave -> 7th -> 6th -> 5th -> 4th) uses the same interval sizes
    # as the ascent from the 4th degree to the octave.
    # The ascent from the 4th degree is: F -> G -> Ab -> Bb -> C'.
    # The intervals are: [1 (F-G), 0.5 (G-Ab), 1 (Ab-Bb), 1 (Bb-C')].
    # So the descending intervals are the same sizes in reverse order of the notes,
    # starting from the last ascending interval.
    # Interval 1 (descent): Octave to 7th (size of Bb to C') -> 1.0
    # Interval 2 (descent): 7th to 6th (size of Ab to Bb) -> 1.0
    # Interval 3 (descent): 6th to 5th (size of G to Ab) -> 0.5
    # Interval 4 (descent): 5th to 4th (size of F to G) -> 1.0
    
    descending_intervals = [
        ascending_intervals[6],  # Interval between 7th and Octave
        ascending_intervals[5],  # Interval between 6th and 7th
        ascending_intervals[4],  # Interval between 5th and 6th
        ascending_intervals[3]   # Interval between 4th and 5th
    ]

    # Step 3: Combine and format the final list of 11 intervals.
    total_intervals = ascending_intervals + descending_intervals

    # Create the final output string in the format {n1,n2,n3,...}
    # Using map(str, ...) to convert each number to a string
    # and ",".join(...) to combine them with commas.
    result_string = "{" + ",".join(map(str, total_intervals)) + "}"
    
    print(result_string)

calculate_maqam_intervals()
<<<_solution_>>>