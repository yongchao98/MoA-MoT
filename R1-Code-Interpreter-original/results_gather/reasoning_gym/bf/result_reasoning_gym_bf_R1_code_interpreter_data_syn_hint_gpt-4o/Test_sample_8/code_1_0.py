# Simulating the Brainf*ck program
tape = [0] * 30000
pointer = 0

# Executing the program
pointer += 1  # Move to the second cell
tape[pointer] = 10  # Set the second cell to 10

# Loop to transfer 10 * 11 to the first cell
while tape[pointer] != 0:
    pointer -= 1  # Move to the first cell
    tape[pointer] += 11
    pointer += 1  # Move back to the second cell
    tape[pointer] -= 1

pointer -= 1  # Move back to the first cell
output = []
output.append(chr(tape[pointer]))  # 'n' (110)
tape[pointer] += 1
output.append(chr(tape[pointer]))  # 'o' (111)
tape[pointer] += 3
output.append(chr(tape[pointer]))  # 'r' (114)
tape[pointer] += 2
output.append(chr(tape[pointer]))  # 't' (116)
tape[pointer] -= 12
output.append(chr(tape[pointer]))  # 'h' (104)

print(''.join(output))