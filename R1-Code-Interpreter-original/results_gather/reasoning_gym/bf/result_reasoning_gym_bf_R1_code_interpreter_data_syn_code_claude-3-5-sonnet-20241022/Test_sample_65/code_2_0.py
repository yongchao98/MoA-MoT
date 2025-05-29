def bf_interpreter(code):
    memory = [0] * 30000
    dp = 0
    output = ""
    
    ip = 0
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
                count = 1
                while count > 0:
                    ip += 1
                    if code[ip] == '[': count += 1
                    if code[ip] == ']': count -= 1
        elif cmd == ']':
            if memory[dp] != 0:
                count = 1
                while count > 0:
                    ip -= 1
                    if code[ip] == '[': count -= 1
                    if code[ip] == ']': count += 1
                ip -= 1
        
        ip += 1
    
    print(output)

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.++.---.-----------.>++++[<++++>-]<+.<"
bf_interpreter(code)