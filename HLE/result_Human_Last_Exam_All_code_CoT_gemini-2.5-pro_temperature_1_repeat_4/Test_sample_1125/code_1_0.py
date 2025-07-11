def solve_last_hat_problem():
    """
    Calculates the number of the 12th person after the numbers of the
    first 11 have been revealed.
    """
    # The numbers on the hats are unique integers from 1 to 12.
    all_possible_numbers = set(range(1, 13))

    # For this demonstration, let's assume the revealed numbers are 1 through 11.
    # In a real scenario, these would be discovered through the pairing strategy.
    revealed_numbers = [1, 5, 12, 3, 7, 10, 2, 8, 4, 9, 6]

    print(f"The 11 revealed numbers are: {sorted(revealed_numbers)}")

    # The total sum of numbers from 1 to 12 is known.
    total_sum = sum(all_possible_numbers)

    # The sum of the numbers that have been revealed.
    sum_of_revealed = sum(revealed_numbers)

    # The 12th person's number is the difference.
    last_persons_number = total_sum - sum_of_revealed

    # The problem asks to show the final equation.
    # We construct the string for the equation.
    revealed_numbers_str = " + ".join(map(str, sorted(revealed_numbers)))
    equation = f"{total_sum} - ({revealed_numbers_str}) = {last_persons_number}"

    print("\nThe 12th person can deduce their number with the following logic:")
    print(equation)
    print(f"\nThe number for the last person is: {last_persons_number}")

solve_last_hat_problem()