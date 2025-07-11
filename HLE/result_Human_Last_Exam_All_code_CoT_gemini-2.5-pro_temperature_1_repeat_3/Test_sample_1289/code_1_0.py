import json

def calculate_maqam_intervals():
    """
    Calculates and prints the musical intervals based on the user's description.
    """
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # Notes: C(1), D(2), E(3), F(4), G(5), Ab(6), B(7), C'(8)
    # Intervals:
    # C-D: Whole tone = 1.0
    # D-E: Whole tone = 1.0
    # E-F: Semitone = 0.5
    # F-G: Whole tone = 1.0
    # G-Ab: Semitone = 0.5
    # Ab-B: Augmented second (3 semitones) = 1.5
    # B-C': Semitone = 0.5
    ascending_intervals = [1.0, 1.0, 0.5, 1.0, 0.5, 1.5, 0.5]

    # Step 2: Define the intervals for the descent.
    # The ascent ends on C'. The descent is 4 notes ending on F.
    # The modified upper register is Jins Nahawand on G: (G-A-Bb-C').
    # So the 4 descending notes sung after C' are Bb, A, G, F.
    # Full note sequence: C, D, E, F, G, Ab, B, C', Bb, A, G, F
    # The intervals for the descent are:
    # C' -> Bb: Whole tone (from Jins Nahawand) = 1.0
    # Bb -> A:  Semitone (from Jins Nahawand) = 0.5
    # A -> G:   Whole tone (from Jins Nahawand) = 1.0
    # G -> F:   Whole tone (from the original lower Jins Ajam) = 1.0
    descending_intervals = [1.0, 0.5, 1.0, 1.0]

    # Step 3: Combine the intervals into a single list.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format and print the output string.
    # The format required is {n1,n2,n3,...}
    # We can create this by converting the list to a string and replacing brackets.
    # Using json.dumps ensures proper formatting for floats without trailing zeros where possible.
    # e.g., 1.0 becomes 1, 0.5 stays 0.5
    
    # Custom formatting to ensure .0 for whole numbers and .5 for halves as per example style
    formatted_intervals = []
    for num in all_intervals:
        if num == int(num):
            formatted_intervals.append(str(int(num)))
        else:
            formatted_intervals.append(str(num))

    output_string = "{" + ",".join(formatted_intervals) + "}"
    
    print("The 11 intervals are:")
    print(output_string)

calculate_maqam_intervals()
<<< {1,1,0.5,1,0.5,1.5,0.5,1,0.5,1,1} >>>