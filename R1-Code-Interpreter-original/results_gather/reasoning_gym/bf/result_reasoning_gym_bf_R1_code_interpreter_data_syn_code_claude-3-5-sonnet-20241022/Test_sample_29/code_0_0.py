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
    for i, cmd in enumerate(code):
        if cmd == '[':
            stack.append(i)
        elif cmd == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start
    
    while ip < len(code):
        cmd = code[ip]
        
        if cmd == '>':
            dp += 1
        elif cmd == '<':
            dp -= 1
        elif cmd == '+':
            memory[dp] = (memory[dp] + 1) % 256
        elif cmd == '-':
            memory[dp] = (memory[dp] - 1) % 256
        elif cmd == '.':
            output += chr(memory[dp])
        elif cmd == '[':
            if memory[dp] == 0:
                ip = brackets[ip]
        elif cmd == ']':
            if memory[dp] != 0:
                ip = brackets[ip]
        
        ip += 1
    
    return output

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.------------.++++++++++++.-----------.+++++++++++++.<"
result = bf_interpreter(code)
print(repr(result))