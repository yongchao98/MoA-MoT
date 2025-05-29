def bf_interpreter(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    # Create loop matching
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
                i = brackets[i]
        elif c == ']':
            if memory[pointer] != 0:
                i = brackets[i]
        i += 1
        
    print(repr(output))

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.-----------.-------.+++++++++++.++++++++.<"
bf_interpreter(code)