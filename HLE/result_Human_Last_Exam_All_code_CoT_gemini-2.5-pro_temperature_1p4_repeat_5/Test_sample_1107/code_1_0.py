def solve():
    """
    This function calculates and prints the final list of files.
    """
    initial_files = ["Makefile", "Opps", "X", "Y", "Z"]
    created_files = ["A", "B", "T"]

    # The file 'OPPS' is assumed to be 'Opps' to match the Makefile.
    # The final set of files is the union of the initial and created files.
    final_files = sorted(list(set(initial_files + created_files)))

    print("The final files in the directory are:")
    for file_name in final_files:
        print(file_name)

solve()