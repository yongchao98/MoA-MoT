def solve_sequence():
    """
    Calculates the next number in the sequence based on a discovered pattern.
    """
    last_number = 2352

    # The rule is n_{k+1} = n_k + 6 * m_k
    # where m_k = d3 * d2 + d1 (d_i are digits of n_k)

    # Deconstruct the last number into its digits
    s_num = str(last_number)
    d3 = int(s_num[0])
    d2 = int(s_num[1])
    d1 = int(s_num[2])
    d0 = int(s_num[3])

    # Calculate the multiplier 'm' based on the rule
    m = d3 * d2 + d1

    # Calculate the difference to be added
    diff = 6 * m

    # Calculate the next number in the sequence
    next_number = last_number + diff

    print("The rule is: next_number = current_number + 6 * (first_digit * second_digit + third_digit)")
    print("For the last number, 2352:")
    print(f"The calculation is: {last_number} + 6 * ({d3} * {d2} + {d1}) = {last_number} + 6 * {m} = {last_number} + {diff} = {next_number}")
    print(f"The next number in the sequence is {next_number}.")

solve_sequence()