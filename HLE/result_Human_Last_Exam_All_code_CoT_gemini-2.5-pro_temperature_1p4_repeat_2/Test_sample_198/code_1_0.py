def solve_for_k():
    """
    Calculates K, the highest number of characters that can be decoded from the Ching.
    """
    # System and memory constants from the problem description
    io_buffer_start_address = 9998

    # Program structure constants
    num_instructions = 108
    opcode_size = 1       # in digits
    register_size = 1     # in digits
    operand_size = 4      # in digits

    # Step 1: Calculate the size of the program's code segment.
    instruction_size = opcode_size + register_size + operand_size
    code_segment_size = num_instructions * instruction_size

    # Step 2: Calculate the size of the available data segment.
    # Data is stored after the code and before the I/O buffers.
    data_start_address = code_segment_size
    data_end_address = io_buffer_start_address - 1
    data_segment_size = data_end_address - data_start_address + 1

    # Step 3: Determine the size of the data block for one encrypted character.
    # The architecture with 4-digit registers and operands suggests data is processed in 4-digit chunks.
    encrypted_char_size = operand_size

    # Step 4: Calculate the highest number of characters (K).
    # This is the floor of the total data size divided by the size of one character's data.
    K = data_segment_size // encrypted_char_size
    
    # Print the detailed calculation as requested.
    print("### Calculation for K ###")
    print(f"1. The size of a single instruction is the sum of its parts: {opcode_size}D (opcode) + {register_size}D (register) + {operand_size}D (operand) = {instruction_size} digits.")
    print(f"2. The total size of the code segment for {num_instructions} instructions is: {num_instructions} * {instruction_size} = {code_segment_size} digits.")
    print(f"3. The code segment occupies memory addresses 0000 to {code_segment_size - 1}.")
    print(f"4. The data segment starts at address {code_segment_size} and ends at {io_buffer_start_address - 1} (before the I/O buffer).")
    print(f"5. The total size of the data segment is: {data_end_address} - {code_segment_size} + 1 = {data_segment_size} digits.")
    print(f"6. Assuming each encrypted character is processed as a {encrypted_char_size}-digit block (based on operand/register size).")
    print("\n--- Final Equation ---")
    print(f"K = floor(Total Data Size / Encrypted Character Size)")
    print(f"K = floor({data_segment_size} / {encrypted_char_size})")
    print(f"K = {K}")
    
    # Also printing the equation with all the initial numbers expanded.
    print("\n--- Expanded Final Equation ---")
    final_eq = f"K = (({io_buffer_start_address} - 1) - ({num_instructions} * ({opcode_size} + {register_size} + {operand_size})) + 1) // {operand_size}"
    print(final_eq)
    print(f"K = (({data_end_address}) - ({num_instructions} * ({instruction_size})) + 1) // {encrypted_char_size}")
    print(f"K = (({data_end_address}) - ({code_segment_size}) + 1) // {encrypted_char_size}")
    print(f"K = ({data_segment_size}) // {encrypted_char_size}")
    print(f"K = {K}")

solve_for_k()