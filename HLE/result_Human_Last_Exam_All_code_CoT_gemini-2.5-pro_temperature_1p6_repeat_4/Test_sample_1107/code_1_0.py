import sys

def solve():
    """
    This script determines and prints the final list of files after running 'make all'.

    The logic follows the execution trace of the make command:
    1. Initial files are: ['Makefile', 'OPPS', 'X', 'Y', 'Z']
    2. 'make all' targets: T, Z, X, Opps.
    3. The rule for 'Opps' ('Opps: T Z') is processed first due to dependency from 'T'.
       - A circular dependency (T -> Opps -> T) is detected and broken by make.
       - 'Z' is up-to-date ('Z' is newer than 'Y').
       - The file 'Opps' does not exist (case-sensitive), so its command 'touch T' is run.
       - New file created: T
    4. The rule for 'X' ('X: Y') is processed next as a dependency for 'T'.
       - 'X' is older than 'Y', so its command 'touch B' is run.
       - New file created: B
    5. The rule for 'T' ('T: Opps X') is now processed.
       - Since its dependencies 'Opps' and 'X' were just rebuilt, 'T' is rebuilt.
       - The command 'touch A' is run.
       - New file created: A
    6. The initial files plus the newly created files make up the final list.
    """
    initial_files = ["Makefile", "OPPS", "X", "Y", "Z"]
    created_files = ["A", "B", "T"]
    
    final_files = initial_files + created_files
    final_files.sort()
    
    print("The final list of files in the directory is:")
    for file_name in final_files:
        print(file_name)

solve()

# The final answer is the sorted list of all files present after the command completes.
# <<<A, B, Makefile, OPPS, T, X, Y, Z>>>