def bf_interpreter(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    # Create bracket matching
    bracket_map = {}
    stack = []
    for i, c in enumerate(code):
        if c == '[':
            stack.append(i)
        elif c == ']':
            start = stack.pop()
            bracket_map[start] = i
            bracket_map[i] = start
    
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
        elif c == '[':
            if memory[pointer] == 0:
                i = bracket_map[i]
        elif c == ']':
            if memory[pointer] != 0:
                i = bracket_map[i]
        i += 1
    return output

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.------------.+++++++++++++.-----.--------."
print(bf_interpreter(code))