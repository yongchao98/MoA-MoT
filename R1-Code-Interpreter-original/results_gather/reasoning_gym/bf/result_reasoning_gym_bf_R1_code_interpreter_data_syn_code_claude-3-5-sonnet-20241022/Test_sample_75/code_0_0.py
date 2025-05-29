def bf_interpreter(code):
    # Initialize memory array with 30000 cells
    memory = [0] * 30000
    # Initialize data pointer and instruction pointer
    dp = 0  # Data pointer
    ip = 0  # Instruction pointer
    output = ""
    
    # Process each instruction
    while ip < len(code):
        command = code[ip]
        
        if command == '>':
            dp += 1
        elif command == '<':
            dp -= 1
        elif command == '+':
            memory[dp] = (memory[dp] + 1) % 256
        elif command == '-':
            memory[dp] = (memory[dp] - 1) % 256
        elif command == '.':
            output += chr(memory[dp])
        elif command == '[':
            if memory[dp] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    ip += 1
                    if code[ip] == '[':
                        bracket_count += 1
                    elif code[ip] == ']':
                        bracket_count -= 1
        elif command == ']':
            if memory[dp] != 0:
                bracket_count = 1
                while bracket_count > 0:
                    ip -= 1
                    if code[ip] == '[':
                        bracket_count -= 1
                    elif code[ip] == ']':
                        bracket_count += 1
                continue
        
        ip += 1
    
    return output

code = "[-]>[-]<>++++++++++[<+++++++++++>-]<+.+++++..---------------.+++++++++++++.<"
result = bf_interpreter(code)
print(repr(result))  # Using repr to show exact output including newlines