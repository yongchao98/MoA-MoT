def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a two-level pattern.
    """
    s = [2, 11, 23, 51, 119]

    print("Analyzing the sequence:", s)
    print("Let's test the pattern a_n = 3 * a_(n-1) +/- constant.")
    print("-" * 30)

    # Establish the pattern for each term
    # Term 2: 11 = 3 * 2 + 5
    print(f"Equation for the 2nd term: {s[1]} = 3 * {s[0]} + 5")
    
    # Term 3: 23 = 3 * 11 - 10
    c1 = 10
    print(f"Equation for the 3rd term: {s[2]} = 3 * {s[1]} - {c1}")

    # Term 4: 51 = 3 * 23 - 18
    c2 = 18
    print(f"Equation for the 4th term: {s[3]} = 3 * {s[2]} - {c2}")

    # Term 5: 119 = 3 * 51 - 34
    c3 = 34
    print(f"Equation for the 5th term: {s[4]} = 3 * {s[3]} - {c3}")
    print("-" * 30)

    # Analyze the sequence of subtracted numbers: c = [10, 18, 34]
    print("The sequence of subtracted numbers is [10, 18, 34, ...]")
    
    # Find the pattern in the differences of c
    diff1 = c2 - c1
    diff2 = c3 - c2
    print(f"The differences in this sequence are {c2}-{c1}={diff1} and {c3}-{c2}={diff2}.")
    print("This is a geometric progression where the difference doubles each time (8, 16, ...).")

    # Predict the next difference
    next_diff = diff2 * 2
    print(f"The next difference will be {diff2} * 2 = {next_diff}.")
    
    # Predict the next number to be subtracted
    c_next = c3 + next_diff
    print(f"The next number to subtract will be {c3} + {next_diff} = {c_next}.")
    print("-" * 30)

    # Calculate the final answer
    last_term = s[-1]
    next_term_in_sequence = 3 * last_term - c_next
    
    print("The final calculation is:")
    print(f"Next Term = 3 * {last_term} - {c_next} = {next_term_in_sequence}")

solve_sequence()