import math

# Step 1: Calculate the total number of distinct instructions available.
# Instructions are: stop (1), turn left (2 levels), turn right (2 levels),
# forward (4 levels), and backward (2 levels).
num_stop_instr = 1
num_left_instr = 2
num_right_instr = 2
num_forward_instr = 4
num_backward_instr = 2

total_distinct_instr = num_stop_instr + num_left_instr + num_right_instr + num_forward_instr + num_backward_instr

print(f"1. Total number of distinct instructions:")
print(f"   {num_stop_instr} (stop) + {num_left_instr} (left) + {num_right_instr} (right) + {num_forward_instr} (forward) + {num_backward_instr} (backward) = {total_distinct_instr} instructions")
print("-" * 20)

# Step 2: Calculate the minimum bits needed to store one instruction.
# Using the formula: bits = ceil(log2(N)), where N is the number of states.
bits_per_instr = math.ceil(math.log2(total_distinct_instr))
print(f"2. Minimum bits per instruction:")
print(f"   ceil(log2({total_distinct_instr})) = {bits_per_instr} bits")
print("-" * 20)

# Step 3: Calculate the total number of instructions needed for one Martian day (sol).
# A sol is 24 hours, 39 minutes, and 35.244 seconds long in Earth time.
# An instruction is issued each Earth minute.
sol_hours = 24
sol_minutes = 39
sol_seconds = 35.244

total_minutes_in_sol = (sol_hours * 60) + sol_minutes + (sol_seconds / 60)
# We need to cover every minute that starts during the sol.
num_instructions_per_day = math.ceil(total_minutes_in_sol)

print(f"3. Total instructions for a Martian Day:")
print(f"   A Martian day is 24h 39m 35.244s, which is {total_minutes_in_sol:.4f} Earth minutes.")
print(f"   Total instructions needed = ceil({total_minutes_in_sol:.4f}) = {num_instructions_per_day} instructions")
print("-" * 20)


# Step 4: Calculate the total memory required in bytes.
# Total bits = total instructions * bits per instruction
# Total bytes = Total bits / 8
total_bits = num_instructions_per_day * bits_per_instr
total_bytes = math.ceil(total_bits / 8)

print(f"4. Total memory size calculation:")
print(f"   Total bits = {num_instructions_per_day} instructions * {bits_per_instr} bits/instruction = {total_bits} bits")
print(f"   Total bytes = {total_bits} bits / 8 bits/byte = {total_bytes} bytes")
print("-" * 20)

print(f"The smallest size of memory required is {total_bytes} bytes.")