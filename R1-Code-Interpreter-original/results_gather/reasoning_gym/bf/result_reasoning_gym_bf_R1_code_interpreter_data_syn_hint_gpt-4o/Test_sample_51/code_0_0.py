# Initialize the memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
tape[pointer] = 0
pointer += 1
tape[pointer] = 0
pointer -= 1
pointer += 1
tape[pointer] += 8
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 15
    pointer += 1
    tape[pointer] -= 1
pointer -= 1
tape[pointer] -= 2
output = chr(tape[pointer])
pointer += 1
tape[pointer] += 3
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] -= 6
    pointer += 1
    tape[pointer] -= 1
output += chr(tape[pointer])
pointer += 1
tape[pointer] += 4
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 5
    pointer += 1
    tape[pointer] -= 1
output += chr(tape[pointer])
pointer += 1
tape[pointer] += 3
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] -= 6
    pointer += 1
    tape[pointer] -= 1
output += chr(tape[pointer])
tape[pointer] += 2
output += chr(tape[pointer])

print(output)