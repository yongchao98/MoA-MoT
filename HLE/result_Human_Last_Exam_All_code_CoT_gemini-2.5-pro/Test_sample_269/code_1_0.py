# The Scheme code involves complex control flow using continuations (call/cc).
# A detailed trace reveals the sequence of numbers printed before termination.
# This Python script simply prints that final sequence.
# Each number in the print statement corresponds to a (write ...) call in the Scheme code.

print(1, 1, 2, 1, 3, 1, 4)