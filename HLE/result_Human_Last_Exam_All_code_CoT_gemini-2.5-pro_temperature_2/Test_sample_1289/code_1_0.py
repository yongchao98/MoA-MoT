def solve_music_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the described performance.
    """

    # Step 1: Define the intervals for the ascending Maqam Zanjaran.
    # Jins Zanjaran on tonic: {0.75, 0.75, 0.5}
    # Jins Hijaz on 4th: {0.5, 1.5, 0.5}
    # Connecting tone to octave: {1.0}
    ascending_intervals = [0.75, 0.75, 0.5, 0.5, 1.5, 0.5, 1.0]

    # Step 2: Define the intervals for the descent using the modified scale.
    # The upper register is now Jins Nahawand on the 4th: {1.0, 0.5, 1.0}
    # The interval to complete the octave in this modified scale is 1.5.
    # The descent from note 8 to 4 traces these upper intervals in reverse pitch order.
    # Interval (8->7): 1.5
    # Interval (7->6): 1.0
    # Interval (6->5): 0.5
    # Interval (5->4): 1.0
    descending_intervals = [1.5, 1.0, 0.5, 1.0]

    # Step 3: Combine the interval lists for the full performance.
    all_intervals = ascending_intervals + descending_intervals

    # Step 4: Format the output as specified.
    # The problem asks to output each number in the final list.
    # The format is a comma-separated string inside curly braces.
    interval_strings = [str(i) for i in all_intervals]
    formatted_output = "{" + ",".join(interval_strings) + "}"

    print(formatted_output)

solve_music_intervals()

print("<<<{0.75,0.75,0.5,0.5,1.5,0.5,1.0,1.5,1.0,0.5,1.0}>>>")