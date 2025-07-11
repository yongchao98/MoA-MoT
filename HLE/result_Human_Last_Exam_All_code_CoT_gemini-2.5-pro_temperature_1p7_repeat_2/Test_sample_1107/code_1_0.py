def get_final_files():
    """
    This script calculates and prints the final list of files after the 'make all' command.
    """
    # Start with the initial set of files in the directory.
    initial_files = {'X', 'Y', 'Z', 'OPPS', 'Makefile'}
    
    # Based on the analysis of the Makefile, the following files are created by 'touch' commands.
    # - 'touch T' is run for the 'Opps' target.
    # - 'touch B' is run for the 'X' target.
    # - 'touch A' is run for the 'T' target.
    # - 'touch C' is NOT run because 'Z' is up-to-date.
    # - 'touch H' is NOT run because 'Q' is not a dependency of 'all'.
    created_files = {'T', 'B', 'A'}
    
    # The final list of files is the union of the initial and created sets.
    final_files = initial_files.union(created_files)
    
    # Sort the files alphabetically for a clean, deterministic output.
    sorted_files = sorted(list(final_files))
    
    print("The files in the directory are:")
    for file_name in sorted_files:
        print(file_name)

get_final_files()