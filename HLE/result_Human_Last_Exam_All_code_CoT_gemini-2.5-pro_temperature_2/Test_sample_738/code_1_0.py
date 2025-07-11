import math

def calculate_memory_for_elon():
    """
    Calculates the smallest size of memory in bytes to store instructions
    for a whole Martian day for the rover Elon.
    """
    # Step 1: Calculate the total number of distinct instructions.
    stop_instr = 1
    left_instr = 2
    right_instr = 2
    forward_instr = 4
    backward_instr = 2

    # The equation for the total number of distinct instructions
    total_distinct_instructions = stop_instr + left_instr + right_instr + forward_instr + backward_instr
    print("Step 1: Calculate total unique instructions")
    print(f"Number of distinct instructions = {stop_instr} (stop) + {left_instr} (left) + {right_instr} (right) + {forward_instr} (forward) + {backward_instr} (backward)")
    print(f"= {total_distinct_instructions} unique instructions")
    print("-" * 30)

    # Step 2: Determine the minimum number of bits needed per instruction.
    bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))
    print("Step 2: Calculate bits needed per instruction")
    print(f"Bits per instruction = ceil(log2({total_distinct_instructions}))")
    print(f"= {bits_per_instruction} bits")
    print("-" * 30)

    # Step 3: Calculate the total number of instructions for a whole Martian day.
    # A Martian day (sol) is approximately 24 hours and 39 minutes long.
    hours_in_sol = 24
    minutes_in_sol_extra = 39
    total_minutes_in_sol = (hours_in_sol * 60) + minutes_in_sol_extra
    print("Step 3: Calculate total instructions for a Martian day")
    print(f"Total instructions needed = ({hours_in_sol} hours * 60 mins/hour) + {minutes_in_sol_extra} minutes")
    print(f"= {total_minutes_in_sol} instructions")
    print("-" * 30)

    # Step 4: Calculate the total memory required in bytes.
    total_bits = total_minutes_in_sol * bits_per_instruction
    total_bytes = math.ceil(total_bits / 8)
    print("Step 4: Calculate total memory in bytes")
    print("Total Memory (Bytes) = ceil( (Total Instructions * Bits Per Instruction) / 8 )")
    print(f"Calculation = ceil( ({total_minutes_in_sol} * {bits_per_instruction}) / 8 )")
    print(f"= ceil( {total_bits} / 8 )")
    print(f"= {total_bytes} bytes\n")

    print(f"The smallest memory size required is {total_bytes} bytes.")


calculate_memory_for_elon()
<<<740>>>