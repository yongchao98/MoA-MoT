def bf_interpreter(code):
    # Pre-process the code to remove any non-BF characters
    code = ''.join(c for c in code if c in '><+-.,[]')
    
    # Initialize memory with just enough cells
    mem = [0] * 100
    ptr = 0
    output = []
    
    # Pre-compute bracket pairs for faster loop processing
    brackets = {}
    stack = []
    for i, c in enumerate(code):
        if c == '[':
            stack.append(i)
        elif c == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start

    i = 0
    while i < len(code):
        c = code[i]
        if c == '>': ptr += 1
        elif c == '<': ptr -= 1
        elif c == '+': mem[ptr] = (mem[ptr] + 1) % 256
        elif c == '-': mem[ptr] = (mem[ptr] - 1) % 256
        elif c == '.': output.append(chr(mem[ptr]))
        elif c == '[' and mem[ptr] == 0: i = brackets[i]
        elif c == ']' and mem[ptr] != 0: i = brackets[i]
        i += 1
    
    print(''.join(output))

code = "[-]>[-]<>++++++++++[<++++++++++>-]<.+++++.+++++++++.++.+++++.<"
bf_interpreter(code)