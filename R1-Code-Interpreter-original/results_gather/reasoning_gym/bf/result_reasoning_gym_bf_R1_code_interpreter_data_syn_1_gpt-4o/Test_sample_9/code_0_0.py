# Initialize memory tape and data pointer
tape = [0] * 30000
ptr = 0

# Execute the Brainf*ck code
tape[ptr] = 0
ptr += 1
tape[ptr] = 0
ptr -= 1
ptr += 1
tape[ptr] += 9
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] += 13
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
tape[ptr] -= 2
output = []
output.append(chr(tape[ptr]))
tape[ptr] -= 6
output.append(chr(tape[ptr]))
tape[ptr] -= 4
output.append(chr(tape[ptr]))
tape[ptr] += 9
output.append(chr(tape[ptr]))
tape[ptr] -= 7
output.append(chr(tape[ptr]))
ptr -= 1

# Print the output
print(''.join(output))