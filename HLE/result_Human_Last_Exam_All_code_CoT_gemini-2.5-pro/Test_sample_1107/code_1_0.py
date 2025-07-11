# This script calculates and prints the final list of files
# in the directory after running 'make all', based on a detailed
# analysis of the Makefile's execution flow.

# 1. Start with the initial set of files present in the directory.
initial_files = {"X", "Y", "Z", "OPPS", "Makefile"}

# 2. Determine which new files are created based on the make process.
# - Rule 'Opps: T Z' executes 'touch T', creating file 'T'.
# - Rule 'X: Y' executes 'touch B', creating file 'B'.
# - Rule 'T: Opps X' executes 'touch A', creating file 'A'.
# - Rule 'Z: Y' is not triggered, so 'C' is not created.
newly_created_files = {"A", "B", "T"}

# 3. Combine the initial and new files to get the final set.
final_files = initial_files.union(newly_created_files)

# 4. Print the final list of all files, sorted alphabetically.
print("The files in the directory after running 'make all' are:")
for file_name in sorted(list(final_files)):
    print(file_name)
