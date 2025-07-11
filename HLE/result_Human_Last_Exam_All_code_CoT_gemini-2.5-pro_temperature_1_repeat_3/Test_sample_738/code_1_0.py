import math

# Step 1: Define the number of options for each instruction type.
stop_instr = 1
turn_left_instr = 2
turn_right_instr = 2
forward_instr = 4
backward_instr = 2

# Sum them to get the total number of distinct instructions.
total_distinct_instructions = stop_instr + turn_left_instr + turn_right_instr + forward_instr + backward_instr

# Step 2: Calculate the minimum bits needed to represent one instruction.
# This is the ceiling of log base 2 of the number of instructions.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

# Step 3: Define the duration of a Martian day in Earth time.
# A Martian day (sol) is 24 hours, 39 minutes, and 35.244 seconds.
hours = 24
minutes = 39
seconds = 35.244

# Calculate the total number of Earth minutes in a Martian day.
total_minutes_in_sol = hours * 60 + minutes + seconds / 60

# An instruction is given each minute, so we need to store one for each minute.
# We round up (take the ceiling) to ensure the entire duration is covered.
total_instructions_for_day = math.ceil(total_minutes_in_sol)

# Step 4: Calculate total bits and convert to bytes.
bits_in_byte = 8
total_memory_bytes = (total_instructions_for_day * bits_per_instruction) / bits_in_byte

# Final output: Print the equation with all the numbers used in the calculation.
print("1. Total distinct instructions:")
print(f"   {stop_instr} (stop) + {turn_left_instr} (left) + {turn_right_instr} (right) + {forward_instr} (forward) + {backward_instr} (backward) = {total_distinct_instructions} instructions")

print("\n2. Bits per instruction:")
print(f"   ceil(log2({total_distinct_instructions})) = {bits_per_instruction} bits")

print("\n3. Total instructions for one Martian day:")
print(f"   ceil({hours}h * 60m/h + {minutes}m + {seconds:.3f}s / 60s/m) = ceil({total_minutes_in_sol:.3f}) = {total_instructions_for_day} instructions")

print("\n4. Final Calculation for memory size in bytes:")
print(f"   ({total_instructions_for_day} total instructions * {bits_per_instruction} bits/instruction) / {bits_in_byte} bits/byte = {int(total_memory_bytes)} bytes")