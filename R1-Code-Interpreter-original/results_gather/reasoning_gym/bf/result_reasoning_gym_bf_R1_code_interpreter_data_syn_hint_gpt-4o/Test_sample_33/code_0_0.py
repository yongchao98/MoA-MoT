# Initialize memory and pointer
memory = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
memory[pointer] = 0
pointer += 1
memory[pointer] = 0
pointer -= 1
pointer += 1
memory[pointer] += 8
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 14
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
memory[pointer] += 2
print(chr(memory[pointer]), end='')  # Output 'r'
memory[pointer] -= 14
print(chr(memory[pointer]), end='')  # Output 'd'
memory[pointer] += 13
print(chr(memory[pointer]), end='')  # Output 'q'
memory[pointer] += 3
print(chr(memory[pointer]), end='')  # Output 't'
memory[pointer] -= 7
print(chr(memory[pointer]), end='')  # Output 'm'
pointer -= 1