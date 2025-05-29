# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
tape[pointer] = 0
pointer += 1
tape[pointer] = 0
pointer += 1
tape[pointer] = 0
pointer -= 1
tape[pointer] += 7

while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 14
    pointer += 1
    tape[pointer] -= 1

pointer -= 1
output = []
output.append(chr(tape[pointer]))
tape[pointer] += 10
output.append(chr(tape[pointer]))
tape[pointer] -= 11
output.append(chr(tape[pointer]))
tape[pointer] += 3
output.append(chr(tape[pointer]))
tape[pointer] += 1
output.append(chr(tape[pointer]))

# Print the output
print(''.join(output))