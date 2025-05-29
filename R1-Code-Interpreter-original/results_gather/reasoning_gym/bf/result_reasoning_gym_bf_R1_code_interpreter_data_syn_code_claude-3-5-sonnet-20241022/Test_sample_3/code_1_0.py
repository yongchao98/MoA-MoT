def bf_interpreter(code):
    tape = [0] * 30000
    pointer = 0
    output = ""
    
    # For debugging
    def print_state():
        return f"Pointer: {pointer}, Current value: {tape[pointer]}, Char: {chr(tape[pointer]) if 0 <= tape[pointer] <= 127 else '?'}"
    
    i = 0
    while i < len(code):
        c = code[i]
        if c == '>':
            pointer += 1
        elif c == '<':
            pointer -= 1
        elif c == '+':
            tape[pointer] = (tape[pointer] + 1) % 256
        elif c == '-':
            tape[pointer] = (tape[pointer] - 1) % 256
        elif c == '.':
            output += chr(tape[pointer])
        elif c == '[':
            if tape[pointer] == 0:
                count = 1
                while count > 0:
                    i += 1
                    if code[i] == '[': count += 1
                    if code[i] == ']': count -= 1
        elif c == ']':
            if tape[pointer] != 0:
                count = 1
                while count > 0:
                    i -= 1
                    if code[i] == '[': count -= 1
                    if code[i] == ']': count += 1
                i -= 1
        i += 1
    
    print(output)

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<.---.+++++++++.-----------.++.<"
bf_interpreter(code)