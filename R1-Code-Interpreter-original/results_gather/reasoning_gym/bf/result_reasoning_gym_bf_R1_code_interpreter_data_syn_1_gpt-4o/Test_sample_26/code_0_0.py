# Initialize memory and pointer
memory = [0] * 30000
pointer = 0
output = []

# Execute the Brainf*ck code
# [-]>[-]<>+++++++++[<+++++++++++++>-]<+.>++++[<---->-]<-.+++++++++++++.+.--------------.<
memory[pointer] = 0  # [-]
pointer += 1
memory[pointer] = 0  # >[-]
pointer -= 1
pointer += 1  # <>
memory[pointer] += 8  # ++++++++
while memory[pointer] != 0:  # [<+++++++++++++>-]
    pointer -= 1
    memory[pointer] += 13
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
output.append(chr(memory[pointer]))  # <+.
pointer += 1
memory[pointer] += 4  # >++++
while memory[pointer] != 0:  # [<---->-]
    pointer -= 1
    memory[pointer] -= 4
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
memory[pointer] -= 1  # <-.
output.append(chr(memory[pointer]))
memory[pointer] += 15  # +++++++++++++++
output.append(chr(memory[pointer]))  # .
memory[pointer] += 1  # +
output.append(chr(memory[pointer]))  # .
memory[pointer] -= 14  # --------------
output.append(chr(memory[pointer]))  # .
pointer -= 1  # <

# Print the output
print(''.join(output))