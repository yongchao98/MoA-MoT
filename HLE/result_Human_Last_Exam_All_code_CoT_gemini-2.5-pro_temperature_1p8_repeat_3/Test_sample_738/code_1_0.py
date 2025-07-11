import math

def calculate_memory_for_elon():
    """
    Calculates the smallest memory size in bytes to store instructions
    for the Elon rover for a whole day on Mars.
    """
    # Step 1: Calculate the total number of unique instructions.
    # stop(1) + turn_left(2) + turn_right(2) + forward(4) + backward(2)
    num_stop_instructions = 1
    num_left_instructions = 2
    num_right_instructions = 2
    num_forward_instructions = 4
    num_backward_instructions = 2

    total_unique_instructions = (num_stop_instructions +
                                 num_left_instructions +
                                 num_right_instructions +
                                 num_forward_instructions +
                                 num_backward_instructions)

    # Step 2: Determine the minimum bits required per instruction.
    # We need to find the smallest 'n' such that 2^n >= total_unique_instructions.
    # This can be calculated with ceiling(log2(total_unique_instructions)).
    bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

    # Step 3: Calculate the total instructions needed for one Martian day.
    # A day on Mars is 24h, 39m, 35s. Instructions are given each Earth minute.
    # We convert the duration of a Martian day to minutes.
    hours_in_sol = 24
    minutes_in_sol = 39
    seconds_in_sol = 35
    total_minutes_in_sol = (hours_in_sol * 60) + minutes_in_sol + (seconds_in_sol / 60)

    # Since an instruction is given "each minute", we need to cover the entire duration.
    # We use ceiling to round up to the next full minute.
    total_instructions_needed = math.ceil(total_minutes_in_sol)

    # Step 4: Calculate the total required memory in bits.
    total_bits = total_instructions_needed * bits_per_instruction

    # Step 5: Convert bits to bytes (1 byte = 8 bits).
    bits_per_byte = 8
    total_bytes = total_bits / bits_per_byte

    print("--- Elon Rover Memory Calculation ---")
    print(f"1. Total unique instructions: {num_stop_instructions} (stop) + {num_left_instructions} (left) + {num_right_instructions} (right) + {num_forward_instructions} (forward) + {num_backward_instructions} (backward) = {total_unique_instructions}")
    print(f"2. Minimum bits per instruction: {bits_per_instruction} bits (since 2^{bits_per_instruction} = {2**bits_per_instruction}, which is >= {total_unique_instructions})")
    print(f"3. Duration of a Martian day: {hours_in_sol}h {minutes_in_sol}m {seconds_in_sol}s, which is approx {total_minutes_in_sol:.2f} Earth minutes.")
    print(f"4. Total instructions for the day: {total_instructions_needed} (one per minute)")
    print(f"5. Total memory required in bits: {total_instructions_needed} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
    print(f"6. Total memory required in bytes: {total_bits} bits / {bits_per_byte} bits/byte = {int(total_bytes)} bytes")
    print("\n--- Final Equation ---")
    print(f"Total instructions ({total_instructions_needed}) * Bits per instruction ({bits_per_instruction}) / Bits per byte ({bits_per_byte})")
    print(f"({total_instructions_needed} * {bits_per_instruction}) / {bits_per_byte} = {int(total_bytes)}")


calculate_memory_for_elon()
<<<740>>>