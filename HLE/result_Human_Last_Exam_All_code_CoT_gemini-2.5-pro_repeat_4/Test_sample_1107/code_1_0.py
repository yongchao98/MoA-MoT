import collections

def solve_make_simulation():
    """
    This script simulates the 'make all' command based on the provided Makefile
    and initial file state to determine the final list of files in the directory.
    """
    
    # Initial files in the directory.
    initial_files = {"X", "Y", "Z", "OPPS", "Makefile"}
    
    # This set will store the names of any new files created.
    created_files = set()
    
    # The timestamps are simplified to represent which file is older.
    # X (older) < Y (newer)
    # Z (newer) > Y (older)
    x_is_older_than_y = True
    z_is_older_than_y = False
    
    # --- Simulation of 'make' execution ---
    
    # 1. 'make' detects the circular dependency: T -> Opps -> T.
    # It breaks the loop, effectively treating the rule for 'Opps' as 'Opps: Z'.
    
    # 2. Rule 'Z: Y' is evaluated.
    # Since Z is newer than Y, the target 'Z' is up-to-date.
    # The command 'touch C' is NOT executed.
    
    # 3. Rule 'Opps: Z' is evaluated.
    # The target file 'Opps' (lowercase 'o') does not exist.
    # Therefore, its command 'touch T' must be executed.
    created_files.add("T")
    
    # 4. Rule 'X: Y' is evaluated.
    # The file 'X' is older than its dependency 'Y'.
    # Therefore, the target 'X' is out-of-date and its command 'touch B' is executed.
    created_files.add("B")
    
    # 5. Rule 'T: Opps X' is evaluated.
    # The target file 'T' now exists (created in step 3).
    # Since 'T' is newer than its file dependency 'X', it is considered up-to-date.
    # The command 'touch A' is NOT executed.
    
    # The final set of files is the union of initial and created files.
    final_files = sorted(list(initial_files.union(created_files)))
    
    print("The final list of files in the directory is:")
    for file_name in final_files:
        print(file_name)

solve_make_simulation()