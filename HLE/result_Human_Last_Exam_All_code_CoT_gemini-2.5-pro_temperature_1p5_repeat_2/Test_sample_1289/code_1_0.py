def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the described performance.
    """
    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # The scale is composed of Jins Ajam on C (1, 1, 0.5) and Jins Hijazkar on G (0.5, 1.5, 0.5),
    # connected by a whole tone interval (1).
    # Intervals: C-D(1), D-E(1), E-F(0.5), F-G(1), G-Ab(0.5), Ab-B(1.5), B-C'(0.5)
    ascending_intervals = [1, 1, 0.5, 1, 0.5, 1.5, 0.5]

    # Step 2: Define the intervals for the descending modified scale.
    # The descent is from C' to F using Jins Nahawand on F.
    # Jins Nahawand on F implies the notes F, G, Ab, Bb. This means B-natural becomes B-flat.
    # The descending path is C' -> Bb -> Ab -> G -> F.
    # Intervals: C'-Bb(1), Bb-Ab(1), Ab-G(0.5), G-F(1)
    descending_intervals = [1, 1, 0.5, 1]

    # Step 3: Combine the ascending and descending intervals.
    # The octave note C' is sung only once, serving as the end of the ascent and the start of the descent.
    # The total number of intervals is 7 (ascent) + 4 (descent) = 11.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output string as requested.
    # Numbers should be formatted to remove trailing '.0' for integers.
    formatted_intervals = []
    for interval in all_intervals:
        if interval == int(interval):
            formatted_intervals.append(str(int(interval)))
        else:
            formatted_intervals.append(str(interval))
    
    # Join the formatted numbers with commas and enclose in curly braces.
    final_output = "{" + ",".join(formatted_intervals) + "}"
    
    print(final_output)

solve_maqam_intervals()
<<<{"1,1,0.5,1,0.5,1.5,0.5,1,1,0.5,1"}>>>