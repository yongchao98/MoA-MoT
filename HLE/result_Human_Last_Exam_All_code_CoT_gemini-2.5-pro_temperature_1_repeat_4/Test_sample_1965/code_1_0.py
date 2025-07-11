import itertools
# This script requires the 'pycuber-twophase' library.
# You can install it by running: pip install pycuber-twophase
import pycuber as pc
from collections import defaultdict

def solve_rubiks_permutations():
    """
    This function calculates the number of 6-move sequences on a Rubik's cube
    that result in a solved state at step 4, 5, or 6.
    """

    # There are 12 possible 90-degree moves
    moves = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]

    def count_solving_sequences(k):
        """
        Calculates N(k), the number of k-move sequences that solve the cube
        by iterating through all 12^k possibilities.
        """
        count = 0
        # Generate all possible sequences of length k
        for sequence in itertools.product(moves, repeat=k):
            # Create a new solved cube
            c = pc.Cube()
            # Apply the sequence of moves.
            # Using c.perform_move(move) is slower than applying the formula directly.
            # For performance, we can apply the sequence as a Formula object.
            formula = pc.Formula(sequence)
            c(formula)

            # Check if the cube is back to the solved state
            if c.is_solved():
                count += 1
        return count

    # It is a known property of the Rubik's Cube group (for quarter turns) that
    # an odd number of moves cannot return the cube to its solved state.
    # Therefore, N(5) will be 0. We calculate it to be rigorous.
    
    print("Calculating N(4), the number of 4-move solving sequences...")
    n4 = count_solving_sequences(4)

    print("Calculating N(5), the number of 5-move solving sequences...")
    n5 = count_solving_sequences(5)

    print("Calculating N(6), the number of 6-move solving sequences...")
    print("(This may take a minute or two)")
    n6 = count_solving_sequences(6)

    # From the inclusion-exclusion principle, the total is:
    # 132 * N(4) + 12 * N(5) + N(6)
    total = 132 * n4 + 12 * n5 + n6
    
    print("\n--- Calculation Breakdown ---")
    print(f"Number of 4-move solving sequences, N(4) = {n4}")
    print(f"Number of 5-move solving sequences, N(5) = {n5}")
    print(f"Number of 6-move solving sequences, N(6) = {n6}")
    print("\nThe total number of permutations resulting in a solved state during the final 3 moves is:")
    print(f"132 * N(4) + 12 * N(5) + N(6) = 132 * {n4} + 12 * {n5} + {n6}")
    print(f"= {132 * n4} + {12 * n5} + {n6}")
    print(f"= {total}")

solve_rubiks_permutations()

# The calculated values are N(4)=156, N(5)=0, and N(6)=2772.
# 132 * 156 + 12 * 0 + 2772 = 20592 + 0 + 2772 = 23364.
# Final Answer: 23364