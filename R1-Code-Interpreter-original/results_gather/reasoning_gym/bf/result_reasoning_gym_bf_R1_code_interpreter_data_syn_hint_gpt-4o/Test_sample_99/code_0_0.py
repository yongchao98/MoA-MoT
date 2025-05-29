# Initialize the memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
tape[pointer] = 0  # [-]
pointer += 1
tape[pointer] = 0  # >[-]
pointer -= 1
pointer += 1
tape[pointer] += 8  # ++++++++
while tape[pointer] != 0:  # [<+++++++++++++>-]
    pointer -= 1
    tape[pointer] += 13
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
print(chr(tape[pointer]), end='')  # <-
tape[pointer] -= 1
tape[pointer] -= 6  # ------
print(chr(tape[pointer]), end='')  # .
tape[pointer] += 5
print(chr(tape[pointer]), end='')  # +
print(chr(tape[pointer]), end='')  # +
tape[pointer] -= 1
print(chr(tape[pointer]), end='')  # -
pointer -= 1