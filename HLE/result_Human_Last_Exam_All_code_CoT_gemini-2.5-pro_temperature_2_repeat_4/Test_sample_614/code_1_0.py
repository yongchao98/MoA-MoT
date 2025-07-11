def solve_riddle():
    """
    This function solves the number riddle based on a logical deduction
    of the given clues and prints the final ordered sequence.
    """
    # The five natural numbers in the puzzle are 1, 2, 3, 4, 5.
    # We will determine their order and place them into a list.

    # Initialize a list to hold the sequence.
    sequence = [0] * 5

    # Clue: "...the fifth, who... lastly follows"
    # This places the number 5 in the last position (index 4).
    sequence[4] = 5

    # Clue: "Number three... protects the last"
    # This places the number 3 in the position before the last (index 3).
    sequence[3] = 3

    # Clue: "Number 1 is the first..."
    # This places the number 1 in the first position (index 0).
    sequence[0] = 1

    # The remaining numbers are 2 and 4 for the remaining positions (index 1 and 2).
    # Clue: "Number 4... always goes behind" number 2.
    # This means 2 comes before 4.
    sequence[1] = 2
    sequence[2] = 4

    # Now, we print the final equation showing the sequence.
    print("The final sequence is:")
    final_equation = f"{sequence[0]}, {sequence[1]}, {sequence[2]}, {sequence[3]}, {sequence[4]}"
    print(final_equation)

solve_riddle()
<<<1, 2, 4, 3, 5>>>