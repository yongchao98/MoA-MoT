# Simulating the Brainf*ck code to find the ASCII values of the output characters

# Initialize the memory tape and pointer
tape = [0] * 30000
pointer = 0
output = []

# Execute the Brainf*ck code
tape[pointer] = 0  # [-]
pointer += 1
tape[pointer] = 0  # >[-]
pointer -= 1
pointer += 1
tape[pointer] = 0  # <>

# +++++++++
tape[pointer] += 9

# [<++++++++++++>-]
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 12
    pointer += 1
    tape[pointer] -= 1

pointer -= 1
output.append(tape[pointer] + 1)  # <+.
tape[pointer] += 1

tape[pointer] -= 4  # ----.
output.append(tape[pointer])

tape[pointer] += 10  # ++++++++++..
output.append(tape[pointer])
output.append(tape[pointer])

tape[pointer] += 6  # ++++++.
output.append(tape[pointer])

# Print the ASCII values of the output characters
print(output)