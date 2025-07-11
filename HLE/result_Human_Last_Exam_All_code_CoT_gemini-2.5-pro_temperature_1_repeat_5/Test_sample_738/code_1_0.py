import math

def calculate_rover_memory():
    """
    Calculates the minimum memory in bytes to store one Martian day of rover instructions.
    """
    # Step 1: Define the number of states for each action to find the total number of distinct instructions.
    stop_instr = 1
    left_instr = 2
    right_instr = 2
    forward_instr = 4
    backward_instr = 2
    
    total_distinct_instructions = stop_instr + left_instr + right_instr + forward_instr + backward_instr

    # Step 2: Determine the minimum number of bits to represent one instruction.
    bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

    # Step 3: Calculate the total number of instructions needed for a whole Martian day.
    # A Martian day (sol) is 24 hours, 39 minutes, and 35.244 seconds long.
    mars_day_hours = 24
    mars_day_minutes = 39
    mars_day_seconds = 35.244
    
    # Convert Martian day duration to total Earth minutes. An instruction is issued each Earth minute.
    mars_day_in_total_seconds = (mars_day_hours * 3600) + (mars_day_minutes * 60) + mars_day_seconds
    total_instructions_to_store = math.ceil(mars_day_in_total_seconds / 60)

    # Step 4: Calculate the total memory size in bits and then convert to bytes.
    bits_in_a_byte = 8
    total_bits = total_instructions_to_store * bits_per_instruction
    total_bytes = total_bits / bits_in_a_byte

    # --- Output the results ---
    print("Step 1: Calculate total distinct instructions")
    print(f"  Instructions: {stop_instr} (stop) + {left_instr} (left) + {right_instr} (right) + {forward_instr} (forward) + {backward_instr} (backward) = {total_distinct_instructions}")
    
    print("\nStep 2: Calculate bits needed per instruction")
    print(f"  Bits per instruction = ceil(log2({total_distinct_instructions})) = {int(bits_per_instruction)} bits")

    print("\nStep 3: Calculate total instructions for a Martian day")
    print(f"  Minutes in a Martian day = ceil(({mars_day_hours}*3600 + {mars_day_minutes}*60 + {mars_day_seconds}) / 60) = {int(total_instructions_to_store)} instructions")

    print("\nStep 4: Calculate total memory in bytes")
    print("  (Bits per instruction * Total instructions) / Bits per byte")
    print("\nFinal Equation:")
    print(f"({int(bits_per_instruction)} * {int(total_instructions_to_store)}) / {bits_in_a_byte} = {int(total_bytes)}")

calculate_rover_memory()