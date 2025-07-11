def solve_for_k():
    """
    Calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # --- Define constants from the problem description ---
    total_symbols_in_ching = 9999
    num_instructions = 108
    
    # Each instruction's size is 1D (opcode) + 1D (register) + 4D (operand).
    instruction_size_in_digits = 1 + 1 + 4
    
    # The architecture's natural data unit size (registers, operands) is 4 digits.
    # We assume this is the size needed to store the information for one character.
    char_data_size_in_digits = 4

    # --- Step 1: Calculate the size of the program's code segment ---
    code_segment_size = num_instructions * instruction_size_in_digits

    # --- Step 2: Calculate the size of the data segment ---
    # The data segment is the total content of the Ching minus the code.
    data_segment_size = total_symbols_in_ching - code_segment_size

    # --- Step 3: Calculate the maximum number of characters (K) ---
    # This is the integer division of the data space by the size of one character's data.
    k = data_segment_size // char_data_size_in_digits

    # --- Print the detailed calculation ---
    print("This script calculates K, the highest number of characters that can be decoded.")
    print("The final calculation is based on the memory capacity for data storage.")
    print("-" * 50)
    print("Final Equation Steps:")
    
    # Print the equation with all numbers substituted
    print(f"K = (Total Symbols in Ching - (Number of Instructions * Size of One Instruction)) // Size of One Character's Data")
    print(f"K = ({total_symbols_in_ching} - ({num_instructions} * {instruction_size_in_digits})) // {char_data_size_in_digits}")
    print(f"K = ({total_symbols_in_ching} - {code_segment_size}) // {char_data_size_in_digits}")
    print(f"K = {data_segment_size} // {char_data_size_in_digits}")
    
    print("-" * 50)
    print(f"The highest number of characters that can be decoded is: {k}")


solve_for_k()

<<<2337>>>