def bf_interpreter(code):
    # Initialize memory with more cells and better pointer tracking
    memory = [0] * 30000
    dp = 0  # data pointer
    ip = 0  # instruction pointer
    output = ""
    
    # Build loop map
    loop_start = []
    loop_map = {}
    
    for i, c in enumerate(code):
        if c == '[':
            loop_start.append(i)
        elif c == ']':
            if loop_start:
                start = loop_start.pop()
                loop_map[start] = i
                loop_map[i] = start
    
    while ip < len(code):
        c = code[ip]
        
        if c == '>':
            dp += 1
        elif c == '<':
            dp -= 1
        elif c == '+':
            memory[dp] = (memory[dp] + 1) % 256
        elif c == '-':
            memory[dp] = (memory[dp] - 1) % 256
        elif c == '.':
            output += chr(memory[dp])
        elif c == '[':
            if memory[dp] == 0:
                ip = loop_map[ip]
        elif c == ']':
            if memory[dp] != 0:
                ip = loop_map[ip]
        
        ip += 1
        
        # Ensure pointer stays within bounds
        dp = dp % 30000
    
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<+.>++++[<---->-]<-.+++++++++++++.++++.>++++[<---->-]<-.<"
print(repr(bf_interpreter(code)))