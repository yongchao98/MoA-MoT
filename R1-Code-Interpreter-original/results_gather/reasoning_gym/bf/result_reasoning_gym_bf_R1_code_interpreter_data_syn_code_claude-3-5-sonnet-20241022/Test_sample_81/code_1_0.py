def bf_interpreter(code):
    # Initialize tape and debug info
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # For debugging
    def print_state(step, command):
        print(f"Step {step}: Command={command}, Pointer={pointer}, Value={tape[pointer]}")
    
    # Find matching brackets
    brackets = {}
    temp_stack = []
    for i, char in enumerate(code):
        if char == '[':
            temp_stack.append(i)
        elif char == ']':
            if temp_stack:
                start = temp_stack.pop()
                brackets[start] = i
                brackets[i] = start
    
    step = 0
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
            print_state(step, command)
        elif command == '[':
            if tape[pointer] == 0:
                code_pointer = brackets[code_pointer]
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = brackets[code_pointer]
        
        code_pointer += 1
        step += 1
    
    print(f"\nFinal output: {output}")

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.+.---------------..-."
bf_interpreter(code)