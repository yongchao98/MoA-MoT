def calculate_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """

    # Step 1: Define the ascending intervals for Maqam Zanjaran.
    # Jins Hijaz (0.5, 1.5, 0.5) on the root, a whole tone (1.0) connector,
    # and Jins Hijaz (0.5, 1.5, 0.5) on the 5th.
    ascending_intervals = [0.5, 1.5, 0.5, 1.0, 0.5, 1.5, 0.5]

    # Step 2: Define the descending intervals for the modified scale.
    # The descent is from the octave to the 4th degree, using a scale
    # with Jins Nahawand on the 4th.
    # The descent path is Octave(C) -> Bb -> Ab -> G -> F(4th).
    # C -> Bb = 1.0 (Whole Tone)
    # Bb -> Ab = 1.0 (Whole Tone)
    # Ab -> G = 0.5 (Semitone)
    # G -> F = 1.0 (Whole Tone)
    descending_intervals = [1.0, 1.0, 0.5, 1.0]

    # Step 3: Combine the lists to get all 11 intervals in order.
    all_intervals = ascending_intervals + descending_intervals

    # Format the list into the required string format: {n1,n2,n3,...}
    # The format specifier '.2f' ensures numbers like 1.0 are printed with decimals.
    formatted_intervals = "{" + ",".join(f"{interval:.2f}".replace(".00", ".0") for interval in all_intervals) + "}"
    
    print("The 11 intervals sung by the musician are:")
    print(formatted_intervals)

calculate_maqam_intervals()

# The final answer in the required format
final_answer = "{" + ",".join(str(i) for i in [0.5, 1.5, 0.5, 1.0, 0.5, 1.5, 0.5, 1.0, 1.0, 0.5, 1.0]) + "}"
print(f"\n<<<{final_answer}>>>")