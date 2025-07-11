import os

def solve():
    """
    Calculates the final list of files after running 'make all'.
    """
    # Initial files in the directory
    initial_files = ['X', 'Y', 'Z', 'OPPS', 'Makefile']

    # According to the makefile analysis:
    # 1. The 'T' and 'Opps' targets have a circular dependency, so 'make' fails to build them.
    #    As a result, 'touch A' and 'touch T' are not executed.
    # 2. Target 'Z' is newer than its dependency 'Y', so 'touch C' is not executed.
    # 3. Target 'X' is older than its dependency 'Y', so 'touch B' is executed.
    # The only new file created is 'B'.
    
    newly_created_file = 'B'

    # Combine the initial files with the newly created one
    final_files = initial_files + [newly_created_file]

    # Sort for consistent output
    final_files.sort()

    print("The files in the directory are:")
    for file_name in final_files:
        print(file_name)

solve()