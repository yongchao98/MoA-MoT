def bf_interpreter(code):
    # Initialize memory array and pointer
    memory = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # Create bracket matching dictionary for loops
    bracket_map = {}
    temp_stack = []
    for i, char in enumerate(code):
        if char == '[':
            temp_stack.append(i)
        elif char == ']':
            if temp_stack:
                start = temp_stack.pop()
                bracket_map[start] = i
                bracket_map[i] = start
    
    while code_pointer < len(code):
        command = code[code_pointer]
        
        if command == '>':
            pointer = (pointer + 1) % 30000
        elif command == '<':
            pointer = (pointer - 1) % 30000
        elif command == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif command == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif command == '.':
            output += chr(memory[pointer])
        elif command == '[':
            if memory[pointer] == 0:
                code_pointer = bracket_map[code_pointer]
        elif command == ']':
            if memory[pointer] != 0:
                code_pointer = bracket_map[code_pointer]
        
        code_pointer += 1
    
    return output

code = "[-]>[-]<>++++++++++[<+++++++++++>-]<+.+++++++.>+++[<------->-]<.>++++[<++++>-]<+.+++++++."
print(repr(bf_interpreter(code)))