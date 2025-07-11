import math

def calculate_memory_for_elon():
    """
    Calculates the smallest size of memory in bytes to store enough
    instruction for a whole day on Mars for Elon.
    """

    # Step 1: Define the number of states for each action.
    stop_states = 1
    turn_left_states = 2
    turn_right_states = 2
    forward_states = 4
    backward_states = 2

    # Calculate the total number of unique instructions.
    total_unique_instructions = stop_states + turn_left_states + turn_right_states + forward_states + backward_states

    # Step 2: Calculate the minimum number of bits to represent all unique instructions.
    # This is the ceiling of log base 2 of the total number of states.
    bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

    # Step 3: Calculate the total number of instructions for one Martian day.
    # A Martian day (sol) is approximately 24 hours and 39 minutes.
    # Instructions are given per Earth minute.
    hours_in_sol = 24
    minutes_in_hour = 60
    extra_minutes_in_sol = 39
    total_minutes_in_sol = (hours_in_sol * minutes_in_hour) + extra_minutes_in_sol

    # Step 4: Calculate the total bits needed for a full day of instructions.
    total_bits_needed = bits_per_instruction * total_minutes_in_sol

    # Step 5: Convert total bits to bytes. There are 8 bits in a byte.
    # We must take the ceiling because memory is allocated in whole bytes.
    bits_in_a_byte = 8
    total_bytes_needed = math.ceil(total_bits_needed / bits_in_a_byte)

    # Print the breakdown of the calculation as requested.
    print(f"The final calculation is based on these numbers:")
    print(f"Total unique instructions = {total_unique_instructions}")
    print(f"Bits required per instruction = {bits_per_instruction}")
    print(f"Total instructions in a Martian Day = {total_minutes_in_sol}")
    print(f"Bits in a byte = {bits_in_a_byte}\n")

    print(f"Final Equation:")
    print(f"Required Memory (bytes) = ceil( (Bits Per Instruction * Total Instructions) / Bits Per Byte )")
    print(f"                            = ceil( ({bits_per_instruction} * {total_minutes_in_sol}) / {bits_in_a_byte} )")
    print(f"                            = ceil( {total_bits_needed} / {bits_in_a_byte} )")
    print(f"                            = {total_bytes_needed}")


if __name__ == '__main__':
    calculate_memory_for_elon()
    # Final answer in the specified format
    # The calculation is ceil( (ceil(log2(1+2+2+4+2))) * (24*60+39) ) / 8 )
    # ceil( (ceil(log2(11))) * (1479) ) / 8 ) = ceil( (4 * 1479) / 8 ) = ceil(5916 / 8) = ceil(739.5) = 740
    print("\n<<<740>>>")
