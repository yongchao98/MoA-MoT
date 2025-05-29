def bf_interpreter(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    # First, let's print the actual values being output
    i = 0
    while i < len(code):
        if code[i] == '[':
            if memory[pointer] == 0:
                count = 1
                while count > 0:
                    i += 1
                    if code[i] == '[': count += 1
                    if code[i] == ']': count -= 1
            else:
                stack = []
                stack.append(i)
        elif code[i] == ']':
            if memory[pointer] != 0:
                i = stack[-1]
                continue
            else:
                stack.pop()
        elif code[i] == '>':
            pointer += 1
        elif code[i] == '<':
            pointer -= 1
        elif code[i] == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif code[i] == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif code[i] == '.':
            output += chr(memory[pointer])
            print(f"Output character: {chr(memory[pointer])} (ASCII: {memory[pointer]})")
        i += 1
    print("\nFinal output:", output)

code = "[-]>[-]<>++++++++[<++++++++++++>-]<+.+++.+++++++++.----.+++++++++++.<"
bf_interpreter(code)