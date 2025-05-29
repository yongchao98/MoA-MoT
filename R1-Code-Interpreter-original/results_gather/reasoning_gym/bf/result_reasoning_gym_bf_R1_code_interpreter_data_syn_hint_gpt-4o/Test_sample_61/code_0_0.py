# Initialize memory and pointer
memory = [0] * 30000
pointer = 0
output = []

# Execute the Brainf*ck code
# >[-]>[-]<>++++++++[<+++++++++++++>-]<.---.----.>+++[<+++++++>-]<.>++++[<---->-]<-.<
pointer += 1
memory[pointer] = 0
pointer += 1
memory[pointer] = 0
pointer -= 1
pointer += 1
memory[pointer] += 8
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 13
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
output.append(chr(memory[pointer]))
memory[pointer] -= 3
output.append(chr(memory[pointer]))
memory[pointer] -= 4
output.append(chr(memory[pointer]))
pointer += 1
memory[pointer] += 3
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 7
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
output.append(chr(memory[pointer]))
pointer += 1
memory[pointer] += 4
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] -= 4
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
memory[pointer] -= 1
output.append(chr(memory[pointer]))
pointer -= 1

# Print the output
print(''.join(output))