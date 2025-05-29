def bf_interpreter(code):
    # Initialize memory array (using a reasonable size)
    memory = [0] * 30000
    pointer = 0
    output = ""
    code_pos = 0
    
    # Helper function to find matching bracket
    def find_matching_bracket(code, start, direction):
        count = 0
        pos = start
        while 0 <= pos < len(code):
            if direction == 1:  # forward
                if code[pos] == '[':
                    count += 1
                elif code[pos] == ']':
                    count -= 1
            else:  # backward
                if code[pos] == '[':
                    count -= 1
                elif code[pos] == ']':
                    count += 1
            if count == 0:
                return pos
            pos += direction
        return -1

    while code_pos < len(code):
        command = code[code_pos]
        
        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif command == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif command == '.':
            output += chr(memory[pointer])
        elif command == '[':
            if memory[pointer] == 0:
                code_pos = find_matching_bracket(code, code_pos, 1)
        elif command == ']':
            if memory[pointer] != 0:
                code_pos = find_matching_bracket(code, code_pos, -1)
        
        code_pos += 1

    print(output)

# Execute the BF code
code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.++.++++..++++++."
bf_interpreter(code)