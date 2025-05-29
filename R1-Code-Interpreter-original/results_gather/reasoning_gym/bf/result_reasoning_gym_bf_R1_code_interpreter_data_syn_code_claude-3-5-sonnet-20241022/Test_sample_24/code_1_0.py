def bf_interpreter(code):
    memory = [0] * 30000
    dp = 0
    ip = 0
    output = ""
    
    # Find matching brackets
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
    
    # Debug: Print memory state at key points
    print(f"Final memory state at relevant cells: {memory[:5]}")
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<-.-----.+++.+.----.<"
result = bf_interpreter(code)
print(f"Character codes: {[ord(c) for c in result]}")
print(f"Result: {result}")