# Initialize memory tape and pointer
tape = [0] * 30000
ptr = 0

# Execute the Brainf*ck code
tape[ptr] = 0
ptr += 1
tape[ptr] = 0
ptr -= 1
ptr += 1
ptr -= 1
tape[ptr] += 9
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] += 13
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
output = chr(tape[ptr])
tape[ptr] += 5
output += chr(tape[ptr])
tape[ptr] -= 7
output += chr(tape[ptr])
ptr += 1
tape[ptr] += 3
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] += 6
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
output += chr(tape[ptr])
tape[ptr] -= 6
output += chr(tape[ptr])

print(output)