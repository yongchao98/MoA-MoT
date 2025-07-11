def solve_maqam_intervals():
    """
    Calculates the sequence of musical intervals based on the problem description.
    """

    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # Jins Hijaz = Semitone (0.5), Augmented Second (1.5), Semitone (0.5)
    # The Maqam is Jins Hijaz + Whole Tone Connector + Jins Hijaz
    jins_hijaz = [0.5, 1.5, 0.5]
    connecting_interval = [1.0]
    ascending_intervals = jins_hijaz + connecting_interval + jins_hijaz
    
    # Step 2: Define the intervals for the descending path with a modified scale.
    # The upper register is modified with Jins Nahawand on the 4th degree.
    # The descent path is: 8 -> 7' -> 6' -> 5' -> 4.
    
    # Jins Nahawand = Whole (1.0), Half (0.5), Whole (1.0).
    # This defines the intervals for notes 4 -> 5' -> 6' -> 7'.
    jins_nahawand_intervals = [1.0, 0.5, 1.0]
    
    # To complete the upper scale, we assume a whole tone from the 7th to the octave (8).
    # So the intervals from note 4 upwards are: I(4,5'), I(5',6'), I(6',7'), I(7',8)
    upper_nahawand_scale_intervals = jins_nahawand_intervals + [1.0]

    # The descending intervals are the same size as the ascending ones, in reverse note order.
    # Interval(8 -> 7') is the size of Interval(7' -> 8), which is the last in the list above.
    # Interval(7' -> 6') is the size of Interval(6' -> 7'), which is the second to last, etc.
    interval_8_to_7prime = upper_nahawand_scale_intervals[3]  # 1.0
    interval_7prime_to_6prime = upper_nahawand_scale_intervals[2] # 1.0
    interval_6prime_to_5prime = upper_nahawand_scale_intervals[1] # 0.5
    interval_5prime_to_4 = upper_nahawand_scale_intervals[0]   # 1.0
    
    descending_intervals = [
        interval_8_to_7prime,
        interval_7prime_to_6prime,
        interval_6prime_to_5prime,
        interval_5prime_to_4
    ]

    # Step 3: Combine both lists and format the output.
    all_intervals = ascending_intervals + descending_intervals

    # Create the final string in the format {x,y,z,...}
    # The map function converts each number in the list to a string.
    # The join method combines them with a comma and space.
    result_string = "{" + ",".join(map(str, all_intervals)) + "}"
    
    print(result_string)

solve_maqam_intervals()
<<<{"0.5,1.5,0.5,1.0,0.5,1.5,0.5,1.0,1.0,0.5,1.0"}>>>