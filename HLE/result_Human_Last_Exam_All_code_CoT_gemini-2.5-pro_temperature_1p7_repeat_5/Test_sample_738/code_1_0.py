import math

# Step 1: Calculate the total number of unique instructions
stop_states = 1
turn_left_states = 2
turn_right_states = 2
forward_states = 4
backward_states = 2

total_unique_instructions = stop_states + turn_left_states + turn_right_states + forward_states + backward_states
print(f"1. Total number of unique instructions:")
print(f"   {stop_states} (stop) + {turn_left_states} (turn left) + {turn_right_states} (turn right) + {forward_states} (forward) + {backward_states} (backward) = {total_unique_instructions} unique instructions")
print("-" * 20)

# Step 2: Determine the number of bits needed to store one instruction
# We need to find the smallest integer 'b' such that 2^b >= total_unique_instructions
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))
print(f"2. Bits needed per instruction:")
print(f"   To represent {total_unique_instructions} states, we need ceil(log2({total_unique_instructions})) = {bits_per_instruction} bits")
print("-" * 20)

# Step 3: Calculate the total number of instructions for one Martian day
# A Martian day (sol) is ~24.6 Earth hours. Instructions are given per Earth minute.
hours_in_sol = 24.6
minutes_in_hour = 60
total_instructions_per_day = int(hours_in_sol * minutes_in_hour)
print(f"3. Total instructions per Martian day:")
print(f"   A Martian day has {hours_in_sol} Earth hours, which is {hours_in_sol} * {minutes_in_hour} = {total_instructions_per_day} minutes/instructions.")
print("-" * 20)

# Step 4: Calculate the total memory size in bytes
bits_per_byte = 8
total_memory_in_bytes = (total_instructions_per_day * bits_per_instruction) / bits_per_byte

print(f"4. Final Calculation for memory in bytes:")
print(f"   ({total_instructions_per_day} instructions * {bits_per_instruction} bits/instruction) / {bits_per_byte} bits/byte = {int(total_memory_in_bytes)} bytes")

<<<738>>>