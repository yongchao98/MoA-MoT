def solve():
    """
    This function determines and prints the final list of files
    after running the 'make all' command based on the provided Makefile.
    """
    # The initial set of files in the directory.
    initial_files = {"X", "Y", "Z", "OPPS", "Makefile"}

    # Files created during the make process.
    # 'touch T' is run for the 'Opps' target.
    # 'touch B' is run for the 'X' target.
    # 'touch A' is run for the 'T' target.
    created_files = {"A", "B", "T"}

    # The final set of files is the union of the initial and created files.
    final_files = initial_files.union(created_files)

    # Sort the files alphabetically for a clean, predictable output.
    sorted_files = sorted(list(final_files))

    # Print the final list of files, one per line.
    for file_name in sorted_files:
        print(file_name)

solve()
<<<A
B
Makefile
OPPS
T
X
Y
Z>>>