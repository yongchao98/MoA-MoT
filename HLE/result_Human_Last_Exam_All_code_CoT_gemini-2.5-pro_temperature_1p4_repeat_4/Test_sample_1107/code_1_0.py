def solve():
    """
    This function determines and prints the list of files in the directory
    after the 'make all' command is executed.
    """
    # The initial set of files in the directory.
    initial_files = ["X", "Y", "Z", "OPPS", "Makefile"]

    # The files that are newly created by the make process.
    # 'touch T' is run for the 'Opps' rule.
    # 'touch B' is run for the 'X' rule.
    # 'touch A' is run for the 'T' rule.
    created_files = ["A", "B", "T"]

    # The final set of files is the combination of the initial files and the created ones.
    final_files = sorted(list(set(initial_files + created_files)))

    print("Files in the directory after 'make all':")
    for file_name in final_files:
        print(file_name)

solve()