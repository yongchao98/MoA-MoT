def bf_interpreter(code):
    cells = [0] * 30000
    cell_ptr = 0
    code_ptr = 0
    output = ""
    
    # Pre-compute bracket pairs
    bracket_map = {}
    temp_stack = []
    for i, char in enumerate(code):
        if char == '[':
            temp_stack.append(i)
        elif char == ']':
            if temp_stack:
                start = temp_stack.pop()
                bracket_map[start] = i
                bracket_map[i] = start

    while code_ptr < len(code):
        command = code[code_ptr]
        
        if command == '>':
            cell_ptr += 1
        elif command == '<':
            cell_ptr -= 1
        elif command == '+':
            cells[cell_ptr] = (cells[cell_ptr] + 1) % 256
        elif command == '-':
            cells[cell_ptr] = (cells[cell_ptr] - 1) % 256
        elif command == '.':
            output += chr(cells[cell_ptr])
        elif command == '[':
            if cells[cell_ptr] == 0:
                code_ptr = bracket_map[code_ptr]
        elif command == ']':
            if cells[cell_ptr] != 0:
                code_ptr = bracket_map[code_ptr]
        
        code_ptr += 1
    
    print(output)

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<-.+++++++++++.---.++++++++.-----------.<"
bf_interpreter(code)