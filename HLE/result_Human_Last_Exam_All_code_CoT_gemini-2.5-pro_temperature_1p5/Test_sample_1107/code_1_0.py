def solve():
    """
    This function determines and prints the final list of files in the directory.
    Based on the Makefile analysis:
    - Initial files: ['X', 'Y', 'Z', 'OPPS', 'Makefile']
    - make all runs.
    - Rule 'X: Y' runs because X is older than Y. Command 'touch B' creates file 'B'.
    - Rule 'Z: Y' does not run because Z is newer than Y. File 'C' is not created.
    - A circular dependency exists between 'T' and 'Opps'. make breaks the loop.
    - Rule 'Opps: T Z' runs because file 'Opps' does not exist. Command 'touch T' creates file 'T'.
    - Rule 'T: Opps X' runs because its dependencies ('Opps' and 'X') were updated. Command 'touch A' creates file 'A'.
    - Rule 'Q: T' is never targeted. File 'H' is not created.
    - The final set of files is the initial set plus the newly created ones.
    """
    initial_files = ["Makefile", "OPPS", "X", "Y", "Z"]
    created_files = ["A", "B", "T"]
    
    final_files = initial_files + created_files
    final_files.sort()
    
    print("The files in the directory are:")
    for file_name in final_files:
        print(file_name)

solve()