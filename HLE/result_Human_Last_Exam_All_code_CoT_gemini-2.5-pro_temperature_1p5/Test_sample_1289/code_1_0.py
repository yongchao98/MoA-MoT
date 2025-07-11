import collections

def solve_music_intervals():
    """
    Calculates and prints the sequence of musical intervals as described.
    """
    # Step 1: Define the 7 ascending intervals of Maqam Zanjaran.
    # Pattern: Whole, 3/4, 3/4, Whole, 3/4, 3/4, Half
    # Values:  1, 0.75, 0.75, 1, 0.75, 0.75, 0.5
    ascending_intervals = [1, 0.75, 0.75, 1, 0.75, 0.75, 0.5]

    # Step 2 & 3: Define the 4 descending intervals.
    # The descent is from the 8th note to the 4th note, using a scale
    # modified with Jins Nahawand on the 4th degree.
    # The modified upper scale intervals (from 4th degree up) are based on Nahawand:
    # 4->5 (Whole), 5->6 (Half), 6->7 (Whole), 7->8 (Whole).
    # The descent path is 8->7, 7->6, 6->5, 5->4.
    # The intervals, recorded as positive values, are:
    # Interval(8->7): Whole tone = 1
    # Interval(7->6): Whole tone = 1
    # Interval(6->5): Half tone = 0.5
    # Interval(5->4): Whole tone = 1
    descending_intervals = [1, 1, 0.5, 1]

    # Step 4: Combine the lists.
    all_intervals = ascending_intervals + descending_intervals

    # Helper function to format numbers (e.g., 1.0 becomes "1")
    def format_number(n):
        if n == int(n):
            return str(int(n))
        else:
            return str(n)

    # Format the list into the required string format: {n1,n2,...}
    # Using a list comprehension to apply the number formatting
    formatted_intervals_str = ",".join([format_number(n) for n in all_intervals])
    final_output = f"{{{formatted_intervals_str}}}"

    # Print the final result.
    print(final_output)

solve_music_intervals()
