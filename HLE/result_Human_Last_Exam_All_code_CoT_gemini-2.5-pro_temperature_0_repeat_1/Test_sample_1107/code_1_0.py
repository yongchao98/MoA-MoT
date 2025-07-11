def solve():
    """
    This function calculates and prints the final list of files
    in the directory after the 'make all' command is executed.
    """
    # Initial files in the directory
    initial_files = ["X", "Y", "Z", "OPPS", "Makefile"]

    # Files created by the make process
    # 'touch B' is run for target X (prerequisite of T)
    # 'touch A' is run for target T
    # 'touch T' is run for target Opps (prerequisite of all)
    created_files = ["A", "B", "T"]

    # The file 'C' is not created because Z is newer than Y.
    # The file 'H' is not created because the target Q is never called.

    # Combine the lists and sort for clean output
    final_files = sorted(initial_files + created_files)

    # Print the final list of files
    for file_name in final_files:
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