def bf_interpreter(code):
    # Initialize tape and debug info
    tape = [0] * 30000
    pointer = 0
    output = ""
    debug = []
    
    # Process each character
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
            debug.append(f"Output ASCII {tape[pointer]} -> '{chr(tape[pointer])}'")
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
    
    # Print debug info
    for d in debug:
        print(d)
    print("\nFinal output:", output)

code = "[-]>[-]<>+++++++[<++++++++++++++>-]<.+++.++.++++++++++++++.-------."
bf_interpreter(code)