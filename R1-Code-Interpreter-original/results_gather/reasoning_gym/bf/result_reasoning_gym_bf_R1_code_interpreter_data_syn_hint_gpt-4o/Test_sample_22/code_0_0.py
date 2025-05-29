# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
tape[pointer] = 0  # [-]
pointer += 1
tape[pointer] = 0  # >[-]
pointer -= 1
pointer += 1  # <>
tape[pointer] += 10  # ++++++++++
while tape[pointer] != 0:  # [<+++++++++++>-]
    pointer -= 1
    tape[pointer] += 13
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
tape[pointer] += 1  # <+
print(chr(tape[pointer]), end='')  # .
tape[pointer] += 1
print(chr(tape[pointer]), end='')  # .
tape[pointer] -= 7  # -------
print(chr(tape[pointer]), end='')  # .
tape[pointer] += 12  # ++++++++++++
print(chr(tape[pointer]), end='')  # .
tape[pointer] -= 8  # --------
print(chr(tape[pointer]), end='')  # .

# The output will be printed by the code