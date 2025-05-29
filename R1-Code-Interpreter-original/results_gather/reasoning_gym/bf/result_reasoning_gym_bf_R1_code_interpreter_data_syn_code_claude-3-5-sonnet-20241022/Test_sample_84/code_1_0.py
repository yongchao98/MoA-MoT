def bf_interpreter_debug(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    values = []  # Store numeric values before conversion to chars
    
    i = 0
    while i < len(code):
        c = code[i]
        if c == '>':
            pointer += 1
        elif c == '<':
            pointer -= 1
        elif c == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif c == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif c == '.':
            output += chr(memory[pointer])
            values.append(memory[pointer])
        elif c == '[':
            if memory[pointer] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif c == ']':
            if memory[pointer] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                i -= 1
        i += 1
    
    print("ASCII values:", values)
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<-.+.---.>++++[<---->-]<.+++++++++++++"
result = bf_interpreter_debug(code)
print("Output:", result)