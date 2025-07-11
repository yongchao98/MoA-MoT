# This script calculates the memory usage for the variables
# in a memory-efficient C interpreter for the X++ language.

# Memory size in bytes for each variable.
# The sizes are chosen based on the minimum required to run the longest possible program.

# 1. Memory for the number of statements 'n'.
# The longest program has n=121 statements.
# An unsigned 8-bit integer (0 to 255) is sufficient.
num_statements_mem = 1

# 2. Memory for the program's variable 'x'.
# The value of x will range from -121 to +121.
# A signed 8-bit integer (-128 to 127) is sufficient.
x_variable_mem = 1

# 3. Memory for the character buffer 'c' used by getchar().
# A standard 'char' type is 1 byte.
char_buffer_mem = 1

# 4. Memory for the loop counter 'i'.
# The loop runs n=121 times, so the counter goes from 0 to 120.
# An unsigned 8-bit integer (0 to 255) is sufficient.
loop_counter_mem = 1

# Calculate the total memory usage for variables.
total_memory = num_statements_mem + x_variable_mem + char_buffer_mem + loop_counter_mem

# Print the final result as a detailed equation.
print("The minimal memory in bytes for the interpreter's variables is the sum of the memory for each variable:")
print(f"Memory for 'n' ({num_statements_mem} byte) + "
      f"Memory for 'x' ({x_variable_mem} byte) + "
      f"Memory for 'c' ({char_buffer_mem} byte) + "
      f"Memory for 'i' ({loop_counter_mem} byte) = {total_memory} bytes")
