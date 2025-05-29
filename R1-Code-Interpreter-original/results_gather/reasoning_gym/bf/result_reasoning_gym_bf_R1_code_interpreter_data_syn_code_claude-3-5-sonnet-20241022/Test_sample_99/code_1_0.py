def bf_interpreter(code):
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # Find matching brackets for loops
    brackets = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start
    
    while code_pointer < len(code):
        command = code[code_pointer]
        
        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            tape[pointer] = (tape[pointer] + 1) % 256
        elif command == '-':
            tape[pointer] = (tape[pointer] - 1) % 256
        elif command == '.':
            output += chr(tape[pointer])
            # Print current value being output
            print(f"Outputting ASCII {tape[pointer]} = '{chr(tape[pointer])}'")
        elif command == '[':
            if tape[pointer] == 0:
                code_pointer = brackets[code_pointer]
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = brackets[code_pointer]
        
        code_pointer += 1
    
    print("\nFinal output:", output)

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<-.------.+++++...-.<"
bf_interpreter(code)