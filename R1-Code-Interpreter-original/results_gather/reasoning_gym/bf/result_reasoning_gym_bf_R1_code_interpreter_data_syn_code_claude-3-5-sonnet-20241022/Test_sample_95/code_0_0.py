def bf_interpreter(code):
    # Initialize memory array with 30000 cells
    memory = [0] * 30000
    # Data pointer
    dp = 0
    # Instruction pointer
    ip = 0
    # Output string
    output = ""
    
    # Loop matching brackets
    bracket_map = {}
    stack = []
    for i, c in enumerate(code):
        if c == '[':
            stack.append(i)
        elif c == ']':
            if stack:
                start = stack.pop()
                bracket_map[start] = i
                bracket_map[i] = start
    
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
                ip = bracket_map[ip]
        elif command == ']':
            if memory[dp] != 0:
                ip = bracket_map[ip]
        
        ip += 1
    
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++>-]<.+++++.-------.++++++++++++.+++."
print(repr(bf_interpreter(code)))