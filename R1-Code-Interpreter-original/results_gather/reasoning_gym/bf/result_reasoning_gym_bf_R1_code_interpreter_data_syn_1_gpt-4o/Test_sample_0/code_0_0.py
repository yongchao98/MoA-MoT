# Initialize memory and pointer
memory = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
code = "[-]>[-]<>++++++++[<+++++++++++++>-]<.---.++++.++++++++++.+.<"
code_pointer = 0
output = []

while code_pointer < len(code):
    command = code[code_pointer]
    
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
                code_pointer += 1
                if code[code_pointer] == '[':
                    open_brackets += 1
                elif code[code_pointer] == ']':
                    open_brackets -= 1
    elif command == ']':
        if memory[pointer] != 0:
            close_brackets = 1
            while close_brackets != 0:
                code_pointer -= 1
                if code[code_pointer] == '[':
                    close_brackets -= 1
                elif code[code_pointer] == ']':
                    close_brackets += 1
    
    code_pointer += 1

print(''.join(output))