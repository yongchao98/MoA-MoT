# Initialize memory and pointer
memory = [0] * 30000
pointer = 0
output = []

# Brainf*ck program
program = "[-]>[-]<>+++++++++[<++++++++++++>-]<.+++++++++.-------.-----------.+++++.<"

# Execute the program
i = 0
while i < len(program):
    command = program[i]
    if command == '>':
        pointer += 1
    elif command == '<':
        pointer -= 1
    elif command == '+':
        memory[pointer] += 1
    elif command == '-':
        memory[pointer] -= 1
    elif command == '.':
        output.append(chr(memory[pointer]))
    elif command == '[':
        if memory[pointer] == 0:
            open_brackets = 1
            while open_brackets != 0:
                i += 1
                if program[i] == '[':
                    open_brackets += 1
                elif program[i] == ']':
                    open_brackets -= 1
    elif command == ']':
        if memory[pointer] != 0:
            close_brackets = 1
            while close_brackets != 0:
                i -= 1
                if program[i] == '[':
                    close_brackets -= 1
                elif program[i] == ']':
                    close_brackets += 1
    i += 1

# Print the output
print(''.join(output))