def solve_wellington_mcq():
    """
    This function analyzes the historical statements about the Duke of Wellington's career
    and identifies the correct options.

    The correct statements identified are:
    1: Relates the Indian mobile commissariat to the successful Peninsular Campaign logistics.
    6: Connects the practice of integrating local forces from India to later British colonial expansion.
    8: Traces the tactical use of flying columns from India to the Peninsula and later to Burma.
    """

    # The numbers of the correct statements, determined by historical analysis.
    correct_options = [1, 6, 8]

    # The problem requires the numbers to be sorted and separated by a comma.
    # The list is already sorted. We will format it for printing.
    answer = ",".join(map(str, correct_options))

    # The final output is printed to the console.
    print(answer)

solve_wellington_mcq()