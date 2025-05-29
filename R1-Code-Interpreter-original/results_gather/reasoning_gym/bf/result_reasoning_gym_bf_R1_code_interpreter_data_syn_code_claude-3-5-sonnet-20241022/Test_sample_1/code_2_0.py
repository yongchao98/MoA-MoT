def bf_debug(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    print("Step by step ASCII values:")
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
            print(f"Output character: '{chr(memory[pointer])}' (ASCII: {memory[pointer]})")
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

code = "[-]>[-]<>++++++++++[<++++++++++>-]<+.-.+++++.---.>+++[<++++++>-]<+.<"
bf_debug(code)