def bf_interpreter(code):
    cells = [0] * 30000  # Create tape with 30000 cells
    cell_pointer = 0
    code_pointer = 0
    output = ""
    
    # Create bracket matching dictionary
    bracket_map = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            if stack:
                start = stack.pop()
                bracket_map[start] = i
                bracket_map[i] = start

    while code_pointer < len(code):
        command = code[code_pointer]
        
        if command == '>':
            cell_pointer += 1
        elif command == '<':
            cell_pointer -= 1
        elif command == '+':
            cells[cell_pointer] = (cells[cell_pointer] + 1) % 256
        elif command == '-':
            cells[cell_pointer] = (cells[cell_pointer] - 1) % 256
        elif command == '.':
            output += chr(cells[cell_pointer])
        elif command == '[':
            if cells[cell_pointer] == 0:
                code_pointer = bracket_map[code_pointer]
        elif command == ']':
            if cells[cell_pointer] != 0:
                code_pointer = bracket_map[code_pointer]
                
        code_pointer += 1
    
    print(output)

# The BF program
program = "[-]>[-]<>+++++++++[<+++++++++++++>-]<+.>++++[<---->-]<-.+++++++++.+++++++.>++++[<---->-]<.<"
bf_interpreter(program)