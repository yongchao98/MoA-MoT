def bf_interpreter(code):
    memory = [0] * 30000
    dp = 0
    output = ""
    
    # For debugging
    def print_relevant_memory():
        return f"Memory at dp={dp}: {memory[max(0,dp-2):dp+3]}"
    
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
            print(f"Output char: {chr(memory[dp])} (ASCII: {memory[dp]})")
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
        
        print(f"After {cmd}: {print_relevant_memory()}")
        ip += 1
    
    return output

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.++.---.-----------.>++++[<++++>-]<+.<"
print(bf_interpreter(code))