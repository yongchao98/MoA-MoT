# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck program
tape[pointer] = 0
pointer += 1
tape[pointer] = 0
pointer -= 1
pointer += 1
tape[pointer] += 9
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 13
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
tape[pointer] -= 2
output = []
output.append(chr(tape[pointer]))
tape[pointer] -= 11
output.append(chr(tape[pointer]))
tape[pointer] -= 7
output.append(chr(tape[pointer]))
tape[pointer] += 11
output.append(chr(tape[pointer]))
tape[pointer] += 8
output.append(chr(tape[pointer]))
pointer -= 1

# Print the output
print(''.join(output))