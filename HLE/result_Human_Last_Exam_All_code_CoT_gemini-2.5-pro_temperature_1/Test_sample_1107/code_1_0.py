# Based on the analysis of the Makefile, we determine the file changes.

# 1. Initial set of files in the directory.
files = {'X', 'Y', 'Z', 'OPPS', 'Makefile'}

# 2. Rule 'X: Y' is triggered because X (10:51) is older than Y (10:52).
#    Command 'touch B' is executed, creating a new file 'B'.
files.add('B')

# 3. Rule 'T: Opps X' is triggered because file 'T' does not exist.
#    Command 'touch A' is executed, creating a new file 'A'.
files.add('A')
#    After the command, `make` also creates the target file 'T'.
files.add('T')

# 4. Rule 'Opps: T Z' is triggered because 'Opps' is now older than the newly created 'T'.
#    The command 'touch T' is executed. This updates T's timestamp but creates no new files
#    as T already exists from the previous step.

# 5. The rule 'Q: T' is not part of the 'all' dependency chain and is therefore ignored.

# The final set of files is the initial set plus the newly created ones.
# We will print them in alphabetical order.
print("The files in the directory after running 'make all' are:")
final_file_list = sorted(list(files))
for file_name in final_file_list:
    print(file_name)