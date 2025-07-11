# This script prints the original, corrected line of code from the shader.
# The bug was a missing semicolon at the end of a variable declaration,
# which is a syntax error that would prevent the shader from compiling.

original_line = "vec4 D = vec4(0.0, 0.5, 1.0, 2.0);// First corner"
print(original_line)