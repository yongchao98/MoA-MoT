def solve_for_k():
    """
    Calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # 1. Define the given constants based on the problem description.
    total_symbols_in_ching = 9999
    num_instructions = 108

    # 2. Determine the size of one instruction in decimal digits (D).
    # Format: [opcode][register][operand] -> 1D + 1D + 4D
    instruction_size_d = 6

    # 3. Determine the size of a single character's representation in decimal digits (D).
    # This is inferred from the machine's architecture (register/operand size).
    char_representation_size_d = 4

    # 4. Calculate the total size of the code section.
    code_section_size = num_instructions * instruction_size_d

    # 5. Calculate the total size of the data section.
    data_section_size = total_symbols_in_ching - code_section_size

    # 6. Calculate K, the maximum number of characters, using integer division.
    K = data_section_size // char_representation_size_d
    
    # 7. Print the final equation with all its numbers and the result, as requested.
    print(f"({total_symbols_in_ching} - ({num_instructions} * {instruction_size_d})) // {char_representation_size_d} = {K}")

solve_for_k()