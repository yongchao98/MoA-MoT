def bf_interpreter(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    # Create bracket matching
    bracket_map = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            start = stack.pop()
            bracket_map[start] = i
            bracket_map[i] = start
    
    i = 0
    while i < len(code):
        command = code[i]
        
        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif command == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif command == '.':
            output += chr(memory[pointer])
        elif command == '[':
            if memory[pointer] == 0:
                i = bracket_map[i]
        elif command == ']':
            if memory[pointer] != 0:
                i = bracket_map[i]
        
        i += 1
    
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++>-]<.+++++.+++++++.------.+++++++++."
print(repr(bf_interpreter(code)))