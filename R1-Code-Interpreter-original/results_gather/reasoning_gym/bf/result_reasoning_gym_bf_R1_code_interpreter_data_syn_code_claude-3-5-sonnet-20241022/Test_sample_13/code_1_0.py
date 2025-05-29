def bf_interpreter(code):
    # Initialize memory and debug info
    memory = [0] * 30000
    pointer = 0
    output = ""
    values = []  # To track ASCII values
    
    # Process each character
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
            values.append(memory[pointer])  # Track ASCII value
        elif c == '[':
            if memory[pointer] == 0:
                count = 1
                while count > 0:
                    i += 1
                    if code[i] == '[': count += 1
                    if code[i] == ']': count -= 1
            else:
                pass
        elif c == ']':
            if memory[pointer] != 0:
                count = 1
                while count > 0:
                    i -= 1
                    if code[i] == '[': count -= 1
                    if code[i] == ']': count += 1
                i -= 1
        i += 1
    
    print(f"Output: {output}")
    print(f"ASCII values: {values}")

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.----.++++++++++..++++++.<"
bf_interpreter(code)