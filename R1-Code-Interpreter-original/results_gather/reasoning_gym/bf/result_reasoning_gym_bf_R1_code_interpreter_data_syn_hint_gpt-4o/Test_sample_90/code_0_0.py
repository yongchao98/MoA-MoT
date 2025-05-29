# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck program
tape[pointer + 1] = 0  # >[-]
tape[pointer + 2] = 0  # >[-]
pointer += 1  # >
pointer -= 1  # <
tape[pointer + 1] += 7  # +++++++
while tape[pointer + 1] != 0:  # [<++++++++++++++>-]
    tape[pointer] += 14
    tape[pointer + 1] -= 1
pointer -= 1  # <
output = chr(tape[pointer])  # .
tape[pointer] += 7  # +++++++
output += chr(tape[pointer])  # .
tape[pointer] -= 5  # -----
output += chr(tape[pointer])  # .
output += chr(tape[pointer])  # .
pointer += 1  # >
tape[pointer] += 3  # +++
while tape[pointer] != 0:  # [<+++++++>-]
    tape[pointer - 1] += 7
    tape[pointer] -= 1
pointer -= 1  # <
output += chr(tape[pointer])  # .
pointer += 1  # <

print(output)