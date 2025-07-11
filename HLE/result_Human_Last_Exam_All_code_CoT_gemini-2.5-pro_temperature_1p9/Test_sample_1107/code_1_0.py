def solve():
    """
    This function determines the final list of files after running 'make all'.
    """
    # Initial files present in the directory
    initial_files = ["X", "Y", "Z", "OPPS", "Makefile"]

    # Files created based on the step-by-step analysis of the Makefile
    # 1. 'touch T' is run because of the 'Opps' target.
    # 2. 'touch B' is run because 'X' is older than its dependency 'Y'.
    # 3. 'touch A' is run because the dependencies for 'T' triggered actions.
    new_files = ["A", "B", "T"]

    # The 'touch C' command for target 'Z' is not run because 'Z' is up-to-date.

    # Combine the initial and new files to get the final list
    final_files = sorted(initial_files + new_files)

    # Print the final list of files
    for file_name in final_files:
        print(file_name)

solve()