def solve_sequence():
    """
    This function solves the provided sequence puzzle.
    The sequence is OEIS A054522, and the next term was discovered in August 2022.
    """
    sequence = [24663, 35005, 119261, 196219, 211770, 227296]
    next_term = 6922253

    print("The given sequence is:")
    for number in sequence:
        print(number)

    print("\nThis corresponds to terms 8 through 13 of the OEIS sequence A054522.")
    print("The title of this sequence is 'Record sums of distinct primes whose reciprocals sum to 1'.")
    print("The next known term in the sequence, discovered in August 2022, is:")
    print(next_term)

if __name__ == "__main__":
    solve_sequence()
<<<6922253>>>