# This code will print the 32-character string that should replace @.
# The original one-liner from the problem description is syntactically flawed,
# but this lambda represents the correct 32-character logic required to solve the problem.
print("lambda s:sorted(set(map(int,s)))")