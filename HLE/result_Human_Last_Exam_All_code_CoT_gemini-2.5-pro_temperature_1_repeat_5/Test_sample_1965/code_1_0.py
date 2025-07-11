import sys
from collections import defaultdict

try:
    import pycuber as pc
except ImportError:
    print("This script requires the pycuber library.")
    print("Please install it using: pip install pycuber")
    sys.exit(1)

def solve_rubiks_problem():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at the 4th, 5th, or 6th move.
    """

    # To make pycuber.Cube objects hashable for use in dictionaries,
    # we can create a simple subclass.
    class HashableCube(pc.Cube):
        def __hash__(self):
            # The __str__ method of a pycuber.Cube provides a canonical
            # representation of its state, which is perfect for hashing.
            return hash(str(self))
        def __eq__(self, other):
            return str(self) == str(other)

    # Define the 12 possible 90-degree moves
    moves = [
        pc.Formula("U"), pc.Formula("U'"), pc.Formula("D"), pc.Formula("D'"),
        pc.Formula("L"), pc.Formula("L'"), pc.Formula("R"), pc.Formula("R'"),
        pc.Formula("F"), pc.Formula("F'"), pc.Formula("B"), pc.Formula("B'"),
    ]

    # The identity state is a solved cube.
    identity_cube = HashableCube()

    # counts is a dictionary mapping a cube state to the number of
    # sequences of moves that produce that state.
    # We start at k=0, where there's 1 way to be in the solved state (by doing nothing).
    counts = {identity_cube: 1}
    
    n_k_I = []
    # We add n_0(I) just to align indices, k=1 to 6.
    n_k_I.append(counts.get(identity_cube, 0))

    # We iterate from k=1 to k=6 to find n_k(I) for each step.
    for k in range(1, 7):
        next_counts = defaultdict(int)
        # For each state reachable in k-1 moves, apply all 12 possible next moves.
        for cube, num_ways in counts.items():
            for move in moves:
                # Create a copy to avoid modifying the cube in the dictionary
                new_cube = cube.copy()
                new_cube.perform_formula(move)
                # Add the number of ways to reach the parent state to the new state's count.
                next_counts[new_cube] += num_ways
        
        counts = next_counts
        n_I = counts.get(identity_cube, 0)
        n_k_I.append(n_I)

    # From the simulation, we get the number of ways to return to identity
    # n_k(I) = 0 for odd k, as expected.
    n_4_I = n_k_I[4]
    n_6_I = n_k_I[6]

    # The final count is given by the formula: 132 * n_4(I) + n_6(I)
    result = 132 * n_4_I + n_6_I

    # Print the final equation and answer
    print(f"The number of ways to return to solved after 4 moves, n_4(I), is: {n_4_I}")
    print(f"The number of ways to return to solved after 6 moves, n_6(I), is: {n_6_I}")
    print("\nThe total number of permutations is calculated by the formula: 132 * n_4(I) + n_6(I)")
    print(f"132 * {n_4_I} + {n_6_I} = {result}")
    
solve_rubiks_problem()