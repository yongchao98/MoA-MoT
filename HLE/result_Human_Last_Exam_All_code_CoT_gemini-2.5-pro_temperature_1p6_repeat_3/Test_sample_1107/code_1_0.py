import collections

# This script simulates the execution of the 'make all' command
# to determine the final list of files in the directory.

# Step 1: Define the initial state of the directory and file timestamps.
# The actual time doesn't matter, only the order.
timestamps = {
    'X': 1,
    'Y': 2,
    'Z': 3,
    'OPPS': 4,
    'Makefile': 5,
}
# The set of files that currently exist.
# Note the case difference: 'OPPS' exists, but the make target is 'Opps'.
existing_files = set(timestamps.keys())

# Step 2: Simulate the 'make all' execution trace.
# The analysis is as follows:
# - Target 'all' depends on 'T', 'Z', 'X', 'Opps'.

# - Evaluating dependency 'X':
#   - Rule 'X: Y'. Since timestamp of X (1) < timestamp of Y (2), X is out of date.
#   - Command 'touch B' is executed.
#   - A new file 'B' is created.
existing_files.add('B')

# - Evaluating dependency 'T':
#   - Rule 'T: Opps X'. 'T' does not exist.
#   - 'make' must first satisfy its dependency 'Opps'.
#   - Rule for 'Opps' is 'Opps: T Z'.
#   - This creates a circular dependency (T -> Opps -> T). 'make' detects this,
#     breaks the loop, and continues.
#   - The file 'Opps' does not exist, so its rule is executed.
#   - Command 'touch T' is executed.
#   - A new file 'T' is created.
existing_files.add('T')
#   - Because 'T' was just created by a dependency's rule, it is considered up-to-date.
#   - The command for 'T' ('touch A') is NOT executed.

# - Evaluating dependency 'Z':
#   - Rule 'Z: Y'. Since timestamp of Z (3) > timestamp of Y (2), Z is up-to-date.
#   - The command 'touch C' is NOT executed.

# - The dependencies 'Opps', 'Q' are processed similarly but don't result
#   in any other new files being created ('Opps' rule runs again but just
#   updates 'T', 'Q' is not a dependency of 'all').

# - Finally, the command for 'all' ('ls') is executed. This command will list
#   all files that exist at this point.

# Step 3: Print the final list of files, sorted alphabetically to mimic 'ls'.
final_file_list = sorted(list(existing_files))

for file_name in final_file_list:
    print(file_name)
