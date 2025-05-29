# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck program
tape[pointer] = 0  # [-]
pointer += 1
tape[pointer] = 0  # >[-]
pointer -= 1
pointer += 1
tape[pointer] += 7  # +++++++
while tape[pointer] != 0:  # [<++++++++++++++>-]
    pointer -= 1
    tape[pointer] += 14
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
output = []
output.append(chr(tape[pointer]))  # <.
pointer += 1
tape[pointer] += 4  # >++++
while tape[pointer] != 0:  # [<++++>-]
    pointer -= 1
    tape[pointer] += 4
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
output.append(chr(tape[pointer]))  # <.
tape[pointer] -= 9  # ---------.
output.append(chr(tape[pointer]))
tape[pointer] -= 6  # ------.
output.append(chr(tape[pointer]))
tape[pointer] += 8  # ++++++++.
output.append(chr(tape[pointer]))

# Print the output
print(''.join(output))