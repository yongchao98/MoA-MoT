def solve_maqam_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """

    # Step 1 & 2: Calculate Ascending Intervals
    # Ascending scale (Maqam Zanjaran as Hijaz on C + Hijaz on G)
    # Pitches in semitones relative to C=0: C=0, Db=1, E=4, F=5, G=7, Ab=8, B=11, C'=12
    ascending_pitches = [0, 1, 4, 5, 7, 8, 11, 12]
    ascending_intervals = []
    for i in range(len(ascending_pitches) - 1):
        # Interval in tones = (semitones / 2)
        interval = (ascending_pitches[i+1] - ascending_pitches[i]) / 2.0
        ascending_intervals.append(interval)

    # Step 3: Calculate Descending Intervals
    # Descending scale (Hijaz on C + Nahawand on F)
    # Notes in descent: C', Bb, Ab, G, F
    # Pitches in semitones: C'=12, Bb=10, Ab=8, G=7, F=5
    descending_path_pitches = [12, 10, 8, 7, 5]
    descending_intervals = []
    for i in range(len(descending_path_pitches) - 1):
        # Interval is always a positive distance
        interval = abs(descending_path_pitches[i+1] - descending_path_pitches[i]) / 2.0
        descending_intervals.append(interval)

    # Step 4: Combine and Format
    all_intervals = ascending_intervals + descending_intervals

    # Format numbers: 1.0 -> '1', 1.5 -> '1.5', etc.
    def format_number(n):
        if n == int(n):
            return str(int(n))
        else:
            return str(n)

    formatted_intervals = [format_number(i) for i in all_intervals]
    
    # Final output string
    result_string = "{" + ",".join(formatted_intervals) + "}"
    print(result_string)

solve_maqam_intervals()
<<<{"0.5,1.5,0.5,1,0.5,1.5,0.5,1,1,0.5,1"}>>>