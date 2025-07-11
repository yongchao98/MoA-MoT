import sys

def solve():
    """
    Solves the grid transformation puzzle.

    The problem describes a 10-step algorithm of rotations, flips, and
    transpositions. By analyzing the geometric transformations, the entire
    sequence can be simplified to a single 90-degree counter-clockwise rotation.

    This script implements that simplified transformation.
    """

    # The initial 11x11 grid provided in the problem
    initial_grid_str = """◫◩▤▤◧◨▥▣▩◨◪
◫◫◫◧◨◪◩▩▩◨▦
▦▧◫▨◧◧◪▥▤▧◫
◧◫▣◩◫◨▨◪▩▤◨
▦▨◪◪▣▧▩▦◨▩▨
◨◫◫◪◪▨▥◪◧▩◨
◧▨▤◩◫◫▣◫◧◨▥
▩▦▥▩◧◧▧▣◪◨◪
◪◨◫▧◫▩▧◧◩▧▩
◩▨▤▨◫▩◨◨◨◫▥
▤▨◫◧◨◪▣▤◨▥◧"""

    # Parse the input string into a list of lists (a 2D array)
    grid = [list(row) for row in initial_grid_str.strip().split('\n')]

    # The entire 10-step algorithm simplifies to one operation:
    # Rotate the grid 90 degrees counter-clockwise.

    # 1. Transpose the grid (swap rows with columns) using zip.
    transposed_grid = zip(*grid)
    
    # 2. Reverse the order of the new rows.
    # We convert each row from a tuple (returned by zip) to a list.
    final_grid = [list(row) for row in transposed_grid][::-1]

    # Print the resulting grid row by row
    for row in final_grid:
        print("".join(row))

solve()