# Initialize memory tape and pointer
tape = [0] * 30000
ptr = 0

# Execute the Brainf*ck code
tape[ptr] = 0
ptr += 1
tape[ptr] = 0
ptr -= 1
ptr += 1
tape[ptr] += 8
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] += 15
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
ptr += 1
tape[ptr] -= 1
output = chr(tape[ptr])
tape[ptr] += 1
output += chr(tape[ptr])
tape[ptr] -= 3
output += chr(tape[ptr])
ptr += 1
tape[ptr] += 4
while tape[ptr] != 0:
    ptr -= 1
    tape[ptr] -= 4
    ptr += 1
    tape[ptr] -= 1
ptr -= 1
ptr += 1
tape[ptr] += 13
# Ensure the value is within the valid range for chr()
if 0 <= tape[ptr] < 256:
    output += chr(tape[ptr])
else:
    output += '?'  # Use '?' for out-of-range values

print(output)