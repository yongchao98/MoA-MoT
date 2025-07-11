import sys

def solve():
    """
    This script calculates and prints the list of files in the directory
    after the 'make all' command is executed, based on the provided Makefile
    and initial file states.
    """
    # Initial files in the directory
    initial_files = ["X", "Y", "Z", "OPPS", "Makefile"]

    # Files created by the make process
    # 'A' is created by the rule for 'T'.
    # 'B' is created by the rule for 'X'.
    # 'T' is created by the rule for 'Opps'.
    # 'C' is NOT created because Z is up-to-date.
    created_files = ["A", "B", "T"]

    # Combine the lists and remove duplicates (though there are none)
    final_files_set = set(initial_files + created_files)

    # Sort the final list of files alphabetically for clean output
    final_files_sorted = sorted(list(final_files_set))

    # Print the final list of files, one per line
    for file_name in final_files_sorted:
        print(file_name)

solve()