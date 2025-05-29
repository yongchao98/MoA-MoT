# Initialize the memory tape and pointer
tape = [0] * 30000
pointer = 0
output = []

# Execute the Brainf*ck code
# [-]>[-]<>++++++++++[<++++++++++>-]<++.-.+++++++++++++++.+.--.<.
tape[pointer] = 0  # [-]
pointer += 1
tape[pointer] = 0  # >[-]
pointer -= 1
pointer += 1  # <>
tape[pointer] += 10  # ++++++++++
while tape[pointer] != 0:  # [<++++++++++>-]
    pointer -= 1
    tape[pointer] += 10
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
tape[pointer] += 2  # <++
output.append(chr(tape[pointer]))  # .
tape[pointer] -= 1  # -
output.append(chr(tape[pointer]))  # .
tape[pointer] += 15  # +++++++++++++++
output.append(chr(tape[pointer]))  # .
tape[pointer] += 1  # +
output.append(chr(tape[pointer]))  # .
tape[pointer] -= 2  # --
output.append(chr(tape[pointer]))  # .
pointer += 1  # <
output.append(chr(tape[pointer]))  # .

# Print the output
print(''.join(output))