# Initialize memory and pointer
memory = [0] * 30000
pointer = 0

# Execute the Brainf*ck program
program = "[-]>[-]<>+++++++++[<++++++++++++>-]<.+++++++++.-------.-----------.+++++.<"
output = []

for command in program:
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
                command = next(program)
                if command == '[':
                    open_brackets += 1
                elif command == ']':
                    open_brackets -= 1
    elif command == ']':
        if memory[pointer] != 0:
            close_brackets = 1
            while close_brackets != 0:
                command = prev(program)
                if command == '[':
                    close_brackets -= 1
                elif command == ']':
                    close_brackets += 1

# Print the output
print(''.join(output))