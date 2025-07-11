def solve_sequence():
    """
    This function solves the number sequence puzzle by finding the underlying pattern.
    """
    seq = [2, 11, 23, 51, 119]

    # Step 1: Find the secondary sequence 'b' based on the pattern x_n+1 = 3*x_n - b_n
    b = []
    for i in range(1, len(seq) - 1):
        b_i = 3 * seq[i] - seq[i+1]
        b.append(b_i)

    # Step 2: Find the pattern in the secondary sequence 'b'
    # The differences in sequence 'b' are 18-10=8, 34-18=16.
    # The difference doubles each time. The next difference is 16 * 2 = 32.
    last_b = b[-1]
    last_diff = b[-1] - b[-2]
    next_diff = last_diff * 2
    next_b = last_b + next_diff

    # Step 3: Calculate the next term in the main sequence
    last_term_in_seq = seq[-1]
    next_term = 3 * last_term_in_seq - next_b

    # Step 4: Print the reasoning and the final equation
    print("The pattern is that each term is 3 times the previous term minus a number.")
    print(f"Let the sequence be x_n. The rule is x_(n+1) = 3 * x_n - b_n for n >= 2.")
    print(f"For x_3 = 23: 3 * {seq[1]} - 10 = 23")
    print(f"For x_4 = 51: 3 * {seq[2]} - 18 = 51")
    print(f"For x_5 = 119: 3 * {seq[3]} - 34 = 119")
    print(f"The sequence of subtracted numbers is {b[0]}, {b[1]}, {b[2]}, ...")
    print(f"The difference in this sequence is {b[1]-b[0]}, {b[2]-b[1]}, ... which doubles each time.")
    print(f"The next difference is {last_diff} * 2 = {next_diff}.")
    print(f"So the next number to subtract is {last_b} + {next_diff} = {next_b}.")
    print("\nFinal calculation:")
    print(f"{3} * {last_term_in_seq} - {next_b} = {next_term}")

solve_sequence()
<<<A>>>