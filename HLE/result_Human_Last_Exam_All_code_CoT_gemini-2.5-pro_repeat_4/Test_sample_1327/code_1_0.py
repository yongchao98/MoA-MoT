def find_next_number():
    """
    This function analyzes the sequence to find the pattern and calculates the next number.
    """
    seq = [2, 11, 23, 51, 119]

    print(f"The given sequence is: {', '.join(map(str, seq))}, ...")
    print("\nLet's analyze the pattern. We find a relationship of the form: next_number = 3 * current_number - s.")

    # Calculate the numbers being subtracted
    s1 = 3 * seq[1] - seq[2]
    s2 = 3 * seq[2] - seq[3]
    s3 = 3 * seq[3] - seq[4]
    
    print(f"1. {seq[2]} = 3 * {seq[1]} - {s1}")
    print(f"2. {seq[3]} = 3 * {seq[2]} - {s2}")
    print(f"3. {seq[4]} = 3 * {seq[3]} - {s3}")

    subtracted_seq = [s1, s2, s3]
    print(f"\nThe sequence of subtracted numbers is: {', '.join(map(str, subtracted_seq))}")

    # Analyze the pattern in the subtracted sequence
    diff1 = s2 - s1
    diff2 = s3 - s2
    print(f"The differences in this new sequence are {diff1} and {diff2}.")
    print("This difference doubles each time (8, 16, ...). The next difference will be 32.")

    # Calculate the next number to be subtracted
    next_s = s3 + (diff2 * 2)
    print(f"The next number to subtract will be {s3} + {diff2 * 2} = {next_s}.")

    # Calculate the final term in the original sequence
    next_term = 3 * seq[4] - next_s
    print("\nTherefore, the next number in the original sequence is calculated as follows:")
    print(f"3 * {seq[4]} - {next_s} = {next_term}")

find_next_number()