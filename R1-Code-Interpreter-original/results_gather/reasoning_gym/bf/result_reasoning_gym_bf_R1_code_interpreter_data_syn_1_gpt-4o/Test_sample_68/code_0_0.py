# Simulating the Brainf*ck code execution
tape = [0] * 30000
ptr = 0

# Execute the Brainf*ck code
tape[ptr] = 0
ptr += 1
tape[ptr] = 0
ptr += 1
tape[ptr] = 0
ptr -= 1
tape[ptr] += 9
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] += 13
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
output = chr(tape[ptr])
ptr += 1
tape[ptr] += 4
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] -= 4
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
tape[ptr] -= 1
output += chr(tape[ptr]) * 2
tape[ptr] += 1
output += chr(tape[ptr])
tape[ptr] += 13
output += chr(tape[ptr])
ptr -= 1

print(output)