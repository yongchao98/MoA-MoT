def solve():
    """
    This script determines the final list of files after running 'make all'.
    """
    # Initial files in the directory
    initial_files = ["Makefile", "OPPS", "X", "Y", "Z"]

    # Files created by the make process based on the logic described.
    # 1. 'touch T' is run for the 'Opps' target.
    # 2. 'touch B' is run for the 'X' target.
    # 3. 'touch A' is run for the 'T' target.
    newly_created_files = ["A", "B", "T"]

    # The file 'C' is never created because the rule for 'Z' is not triggered.
    # The file 'H' is never created because the target 'Q' is not part of 'all'.

    # Combine the lists and sort for a clean, deterministic output.
    final_files = sorted(initial_files + newly_created_files)

    # Print the final list of files, one per line.
    for f in final_files:
        print(f)

solve()