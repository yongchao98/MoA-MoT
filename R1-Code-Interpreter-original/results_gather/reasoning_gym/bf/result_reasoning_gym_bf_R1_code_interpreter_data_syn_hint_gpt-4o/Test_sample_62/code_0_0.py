# Simulating the Brainf*ck program
tape = [0] * 30000
pointer = 0

# Program execution
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
output = chr(tape[pointer])
tape[pointer] -= 3
output += chr(tape[pointer])
tape[pointer] -= 11
output += chr(tape[pointer])
tape[pointer] -= 4
output += chr(tape[pointer])
tape[pointer] += 10
output += chr(tape[pointer])

print(output)