def bf_interpreter(code):
    memory = [0] * 30000
    dp = 0
    ip = 0
    output = ""
    
    # Match brackets
    brackets = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start
    
    # Debug array to track values
    debug_values = []
    
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
            debug_values.append(memory[dp])
        elif command == '[':
            if memory[dp] == 0:
                ip = brackets[ip]
        elif command == ']':
            if memory[dp] != 0:
                ip = brackets[ip]
        
        ip += 1
    
    print(f"ASCII values at output: {debug_values}")
    return output

code = "[-]>[-]<>++++++++[<++++++++++++>-]<+.>++++[<+++++>-]<.>++++[<---->-]<-.+++++.++++++.<"
result = bf_interpreter(code)
print(result)