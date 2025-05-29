# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
# [-]
while tape[pointer] != 0:
    tape[pointer] -= 1

# >[-]
pointer += 1
while tape[pointer] != 0:
    tape[pointer] -= 1

# <>
pointer -= 1

# +++++++
tape[pointer] += 7

# [<+++++++++++++++>-]
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 15
    pointer += 1
    tape[pointer] -= 1

# <.
pointer -= 1
output = chr(tape[pointer])

# -----
tape[pointer] -= 5
output += chr(tape[pointer])

# +++++
tape[pointer] += 5
output += chr(tape[pointer])

# ++++++
tape[pointer] += 6
output += chr(tape[pointer])

# --.<
tape[pointer] -= 2
output += chr(tape[pointer])
pointer -= 1

# Print the output
print(output)