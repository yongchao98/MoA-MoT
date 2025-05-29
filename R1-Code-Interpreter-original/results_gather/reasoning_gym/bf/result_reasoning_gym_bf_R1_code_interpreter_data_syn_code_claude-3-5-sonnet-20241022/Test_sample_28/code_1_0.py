def bf_interpreter_debug(code):
    memory = [0] * 30000
    dp = 0
    ip = 0
    output = ""
    debug_info = []
    
    jump_table = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            if stack:
                start = stack.pop()
                jump_table[start] = i
                jump_table[i] = start
    
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
            debug_info.append(f"Printed ASCII {memory[dp]} -> '{chr(memory[dp])}'")
        elif command == '[':
            if memory[dp] == 0:
                ip = jump_table[ip]
        elif command == ']':
            if memory[dp] != 0:
                ip = jump_table[ip]
        
        ip += 1
    
    print("Debug information:")
    for info in debug_info:
        print(info)
    print("\nFinal output:", repr(output))

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.++++.>+++[<------->-]<-.++++++++++++.+++.<"
bf_interpreter_debug(code)