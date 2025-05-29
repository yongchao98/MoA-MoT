# Initialize memory and pointer
memory = [0] * 30000
pointer = 0
output = []

# Execute the Brainf*ck code
code = "[-]>[-]<>++++++++++[<+++++++++++>-]<+.+++++..---------------.+++++++++++++.<"
i = 0
while i < len(code):
    command = code[i]
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
                if code[i] == '[':
                    open_brackets += 1
                elif code[i] == ']':
                    open_brackets -= 1
    elif command == ']':
        if memory[pointer] != 0:
            close_brackets = 1
            while close_brackets != 0:
                i -= 1
                if code[i] == '[':
                    close_brackets -= 1
                elif code[i] == ']':
                    close_brackets += 1
    i += 1

# Print the output
print(''.join(output))