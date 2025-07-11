import sys

# Based on the analysis of the Makefile, this is the final list of files.
# - Initial files: X, Y, Z, OPPS, Makefile
# - The 'make' command 'touch B' is executed because X is older than Y.
# - Other 'touch' commands are not executed due to either a circular dependency or the target being up-to-date.
# - Therefore, the new file 'B' is created.

final_files = ["X", "Y", "Z", "OPPS", "Makefile", "B"]

# Sort the list for a clean, predictable output.
final_files.sort()

# Print each file name.
for file_name in final_files:
    print(file_name)
