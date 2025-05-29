def bf_interpreter(code):
    # Initialize memory with 30000 cells
    memory = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # Create a dictionary to match brackets for loops
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
    
    while code_pointer < len(code):
        command = code[code_pointer]
        
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
                code_pointer = brackets[code_pointer]
        elif command == ']':
            if memory[pointer] != 0:
                code_pointer = brackets[code_pointer]
        
        code_pointer += 1
    
    return output

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.++.+++++++.-------------.----.<"
print(repr(bf_interpreter(code)))