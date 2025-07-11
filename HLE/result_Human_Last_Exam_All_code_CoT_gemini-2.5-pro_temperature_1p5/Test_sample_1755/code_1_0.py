# The line number with the error.
line_number = 5

# The number of edit operations required for the fix.
# Changing the type 'char' to 'int' is a single operation.
edit_operations = 1

# The problem asks for the answer in the format l:n.
# The following print statement constructs this "final equation"
# by outputting each number separated by a colon.
print(f"{line_number}:{edit_operations}")