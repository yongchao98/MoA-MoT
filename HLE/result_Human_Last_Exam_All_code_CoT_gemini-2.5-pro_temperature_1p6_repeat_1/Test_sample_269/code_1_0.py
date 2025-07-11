# The Scheme code uses continuations to create a complex control flow.
# Tracing the execution reveals the following sequence of `write` operations.
# The code below prints the numbers from that sequence.
# Each number is passed as a separate argument to the print function,
# which corresponds to the separate `(write ...)` calls in the original code.

print(1, 1, 2, 1, 3, 1, 4)