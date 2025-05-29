def bf_interpreter(code):
    tape = [0] * 30000
    ptr = 0
    output = ""
    i = 0
    
    while i < len(code):
        c = code[i]
        if c == '>':
            ptr += 1
        elif c == '<':
            ptr -= 1
        elif c == '+':
            tape[ptr] = (tape[ptr] + 1) % 256
        elif c == '-':
            tape[ptr] = (tape[ptr] - 1) % 256
        elif c == '.':
            output += chr(tape[ptr])
        elif c == '[':
            if tape[ptr] == 0:
                nesting = 1
                while nesting > 0:
                    i += 1
                    if code[i] == '[': nesting += 1
                    if code[i] == ']': nesting -= 1
            else:
                continue
        elif c == ']':
            if tape[ptr] != 0:
                nesting = 1
                while nesting > 0:
                    i -= 1
                    if code[i] == '[': nesting -= 1
                    if code[i] == ']': nesting += 1
                continue
        i += 1
    print(output)

# Let's break down the code:
# [-]>[-]<>++++++++++[<++++++++++>-]<.+++++.+++++++++.++.+++++.<
code = "[-]>[-]<>++++++++++[<++++++++++>-]<.+++++.+++++++++.++.+++++.<"
bf_interpreter(code)