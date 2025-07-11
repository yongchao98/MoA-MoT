# This script simulates the `make all` command to determine the final files
# in the directory.

# 1. Initial files in the directory.
initial_files = ["X", "Y", "Z", "OPPS", "Makefile"]

# 2. Files to be created based on the Makefile logic.
# An empty list to start, we will add files if their creation rules are met.
new_files = []

# 3. Simulate the 'make' execution process.

# Rule for 'X': 'X: Y'. File X (10:51) is older than Y (10:52).
# This means the target is out-of-date and the command must be run.
# The command is `touch B`, so file 'B' is created.
new_files.append("B")

# Rule for 'Z': 'Z: Y'. File Z (10:54) is newer than Y (10:52).
# This means the target is up-to-date. The command `touch C` is not run.

# Rules for 'T' and 'Opps':
# 'T' depends on 'Opps', and 'Opps' depends on 'T'.
# This is a circular dependency. 'make' will detect this error and will not
# execute the commands for either target (`touch A` or `touch T`).

# Rule for 'all':
# Since its dependencies 'T' and 'Opps' cannot be built due to the error,
# the `make` process for 'all' fails, and its command `ls` is not run.

# 4. Determine the final list of files.
final_files = sorted(initial_files + new_files)

# 5. Print the final list of files.
for file_name in final_files:
    print(file_name)
