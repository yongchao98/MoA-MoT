import os

def solve():
    """
    This script determines and prints the final list of files
    after running the 'make all' command as described.
    """
    # Initial files present in the directory
    initial_files = ["Makefile", "OPPS", "X", "Y", "Z"]

    # Files created by the make process based on the trace
    # 'touch T' is run for the 'Opps' target.
    # 'touch B' is run for the 'X' target.
    # 'touch A' is run for the 'T' target.
    # 'touch C' is NOT run because target Z is up-to-date.
    created_files = ["A", "B", "T"]

    # The final set of files is the union of the initial and created files.
    final_files = sorted(list(set(initial_files + created_files)))

    # Print the final list of files, one per line.
    for f in final_files:
        print(f)

solve()