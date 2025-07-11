import sys

def solve():
    """
    This script calculates and prints the list of files in the directory
    after running the 'make all' command as described.
    """

    # Initial files in the directory
    initial_files = ["Makefile", "OPPS", "X", "Y", "Z"]

    # Files created by the make process
    # 'touch B' is run because Y is newer than X.
    # 'touch T' is run because the target 'Opps' does not exist.
    # 'touch A' is run because the target 'T' does not exist initially.
    # 'touch C' is not run because Z is newer than Y.
    created_files = ["A", "B", "T"]

    # The final list of files is the union of the initial and created files.
    final_files = sorted(list(set(initial_files + created_files)))

    print("The final list of files in the directory is:")
    for file_name in final_files:
        print(file_name)

solve()
