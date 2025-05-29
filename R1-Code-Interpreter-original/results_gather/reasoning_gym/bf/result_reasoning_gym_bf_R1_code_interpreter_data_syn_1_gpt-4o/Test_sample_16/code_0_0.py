# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
tape[pointer] = 0
pointer += 1
tape[pointer] = 0
pointer -= 1
pointer += 1
tape[pointer] += 9
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 12
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
tape[pointer] += 1
output = chr(tape[pointer])
tape[pointer] += 1
output += chr(tape[pointer])
tape[pointer] += 4
output += chr(tape[pointer])
output += chr(tape[pointer])
tape[pointer] += 6
output += chr(tape[pointer])
pointer -= 1

# Print the output
print(output)