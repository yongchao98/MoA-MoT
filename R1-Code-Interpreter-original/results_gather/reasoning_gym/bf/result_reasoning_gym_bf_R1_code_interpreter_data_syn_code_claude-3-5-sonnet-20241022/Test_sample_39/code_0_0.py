def bf_interpreter(code):
    # Initialize memory array with 30000 cells
    memory = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    while code_pointer < len(code):
        command = code[code_pointer]
        
        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif command == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif command == '.':
            output += chr(memory[pointer])
        elif command == '[':
            if memory[pointer] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer += 1
                    if code[code_pointer] == '[':
                        bracket_count += 1
                    elif code[code_pointer] == ']':
                        bracket_count -= 1
        elif command == ']':
            if memory[pointer] != 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer -= 1
                    if code[code_pointer] == '[':
                        bracket_count -= 1
                    elif code[code_pointer] == ']':
                        bracket_count += 1
                continue
        
        code_pointer += 1
    
    return output

code = "[-]>[-]<>++++++++++[<++++++++++>-]<.++++++++++++++.---..---.<"
result = bf_interpreter(code)
print(result)