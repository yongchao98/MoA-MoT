def bf_interpreter(code):
    # Initialize memory array with 30000 cells
    memory = [0] * 30000
    # Data pointer
    dp = 0
    # Instruction pointer
    ip = 0
    # Output string
    output = ""
    
    # Loop brackets matching
    brackets = {}
    stack = []
    
    # Find matching brackets first
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start
    
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
                ip = brackets[ip]
        elif command == ']':
            if memory[dp] != 0:
                ip = brackets[ip]
        
        ip += 1
    
    return output

code = "[-]>[-]<>+++++++[<+++++++++++++++>-]<.+++++.-------.++++++++.+++++."
print(repr(bf_interpreter(code)))