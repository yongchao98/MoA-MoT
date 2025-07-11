import math

# Step 1: Define the number of options for each instruction.
stop_levels = 1
turn_left_levels = 2
turn_right_levels = 2
forward_levels = 4
backward_levels = 2

# Calculate the total number of unique instructions.
total_instruction_types = stop_levels + turn_left_levels + turn_right_levels + forward_levels + backward_levels

# Step 2: Calculate the bits needed to represent one instruction.
bits_per_instruction = math.ceil(math.log2(total_instruction_types))

# Step 3: Calculate the total number of instructions needed for a whole Martian day.
# A Martian day is 24 hours, 39 minutes, and 35 seconds.
# An instruction is sent every Earth minute.
# We round up to the next minute because of the extra 35 seconds.
hours_in_sol = 24
minutes_in_sol = 39
total_instructions_per_day = (hours_in_sol * 60) + minutes_in_sol + 1

# Step 4: Calculate the total bits and convert to bytes.
total_bits = total_instructions_per_day * bits_per_instruction
total_bytes = total_bits / 8

# Print the breakdown of the calculation and the final equation.
print(f"1. Total unique instructions = {stop_levels} + {turn_left_levels} + {turn_right_levels} + {forward_levels} + {backward_levels} = {total_instruction_types}")
print(f"2. Bits per instruction = ceil(log2({total_instruction_types})) = {bits_per_instruction}")
print(f"3. Instructions per day = ({hours_in_sol} * 60) + {minutes_in_sol} + 1 = {total_instructions_per_day}")
print("\nFinal calculation:")
print(f"Memory in Bytes = (({hours_in_sol} * 60) + {minutes_in_sol} + 1) * {bits_per_instruction} / 8 = {int(total_bytes)}")