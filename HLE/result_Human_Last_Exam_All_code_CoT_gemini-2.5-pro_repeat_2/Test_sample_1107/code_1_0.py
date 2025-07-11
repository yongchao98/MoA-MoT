import sys

# This script calculates and prints the final list of files after 'make all' is run.

# Step 1: Define the initial set of files in the directory.
initial_files = {"X", "Y", "Z", "OPPS", "Makefile"}

# Step 2: Analyze the Makefile execution to determine which new files are created.
# - The target 'X' is rebuilt because it is older than its dependency 'Y'.
#   Its recipe 'touch B' creates the file 'B'.
# - Because 'X' was rebuilt and 'T' does not exist, the recipe for 'T' ('touch A') is run.
#   This creates the file 'A'.
# - The target 'Opps' is always out-of-date. Its dependency 'T' was just made.
#   Its recipe 'touch T' is run, creating the file 'T'.
newly_created_files = {"A", "B", "T"}

# Step 3: Combine the initial and new files to get the final list.
final_files = sorted(list(initial_files.union(newly_created_files)))

# Step 4: Print the final list of files, one per line.
print("The final files in the directory are:")
for file_name in final_files:
    print(file_name)

# This is a string representation of the final answer for the user.
final_answer_string = ", ".join(final_files)
# Use a specific print format for the final answer to be captured.
sys.stdout.write(f"\n<<<>>>\n{final_answer_string}\n")