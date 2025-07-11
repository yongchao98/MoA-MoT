import sys

def solve():
    """
    This function determines and prints the list of files present in the directory
    after the `make all` command is executed, based on the provided Makefile and
    initial file list.
    """

    # Initial files in the directory.
    # The file is 'OPPS' (uppercase), which is different from the target 'Opps'.
    initial_files = {"X", "Y", "Z", "OPPS", "Makefile"}

    # Files are created based on the step-by-step analysis:
    # - 'A' is created by the rule for target 'T'.
    # - 'B' is created by the rule for target 'X'.
    # - 'T' is created by the rule for target 'Opps'.
    # - 'C' is never created because target 'Z' is always up-to-date.
    created_files = {"A", "B", "T"}

    # The final set of files is the union of the initial and created files.
    final_files = sorted(list(initial_files.union(created_files)))

    # Print the final list of files, one per line.
    for file_name in final_files:
        print(file_name)

solve()