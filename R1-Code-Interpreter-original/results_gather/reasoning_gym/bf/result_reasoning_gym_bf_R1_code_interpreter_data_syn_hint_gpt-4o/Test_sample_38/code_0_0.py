# Simulating the Brainf*ck code execution
memory = [0] * 30000
pointer = 0
output = []

# Initial setup
memory[pointer] = 0
pointer += 1
memory[pointer] = 0
pointer -= 1
pointer += 1
memory[pointer] += 8

# Loop to set the first cell
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 12
    pointer += 1
    memory[pointer] -= 1

# Output sequence
pointer -= 1
memory[pointer] += 1
output.append(chr(memory[pointer]))  # 'a'
memory[pointer] += 3
output.append(chr(memory[pointer]))  # 'd'
memory[pointer] += 9
output.append(chr(memory[pointer]))  # 'm'
memory[pointer] -= 4
output.append(chr(memory[pointer]))  # 'i'
memory[pointer] += 11
output.append(chr(memory[pointer]))  # 'v'

# Print the output
print(''.join(output))