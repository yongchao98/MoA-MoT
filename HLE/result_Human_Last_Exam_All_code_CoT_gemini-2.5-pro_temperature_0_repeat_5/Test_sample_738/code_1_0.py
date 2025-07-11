import math

# Step 1: Calculate the total number of unique instructions.
# Stop: 1, Turn Left: 2, Turn Right: 2, Forward: 4, Backward: 2
num_stop_instr = 1
num_left_instr = 2
num_right_instr = 2
num_forward_instr = 4
num_backward_instr = 2

total_unique_instructions = num_stop_instr + num_left_instr + num_right_instr + num_forward_instr + num_backward_instr
print(f"1. Calculating the number of unique instructions:")
print(f"   Stop (1) + Turn Left ({num_left_instr}) + Turn Right ({num_right_instr}) + Forward ({num_forward_instr}) + Backward ({num_backward_instr}) = {total_unique_instructions} unique instructions.")
print("-" * 20)

# Step 2: Determine the minimum number of bits per instruction.
# We need to find the smallest integer 'n' such that 2^n >= total_unique_instructions.
# This is equivalent to ceil(log2(total_unique_instructions)).
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))
print(f"2. Calculating the bits needed per instruction:")
print(f"   To represent {total_unique_instructions} states, we need ceil(log2({total_unique_instructions})) = {bits_per_instruction} bits per instruction.")
print("-" * 20)

# Step 3: Calculate the total number of instructions for a Martian day.
# A Martian day (sol) is 24 hours, 39 minutes, and 35.244 seconds.
# We convert this to Earth minutes.
hours_in_sol = 24
minutes_in_sol = 39
seconds_in_sol = 35.244
total_minutes_in_sol = (hours_in_sol * 60) + minutes_in_sol + (seconds_in_sol / 60)

# Since an instruction is sent each minute, we need to cover the whole duration.
# We take the ceiling to ensure we have an instruction for the final partial minute.
total_instructions_for_a_day = math.ceil(total_minutes_in_sol)
print(f"3. Calculating the total instructions for a Martian day:")
print(f"   A Martian day is ~{total_minutes_in_sol:.2f} Earth minutes long.")
print(f"   Total instructions needed = ceil({total_minutes_in_sol:.2f}) = {total_instructions_for_a_day} instructions.")
print("-" * 20)

# Step 4: Calculate the total memory size in bits.
total_bits = total_instructions_for_a_day * bits_per_instruction
print(f"4. Calculating the total memory size in bits:")
print(f"   Total bits = {total_instructions_for_a_day} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits.")
print("-" * 20)

# Step 5: Convert the total memory size to bytes.
# 1 byte = 8 bits.
total_bytes = total_bits / 8
print(f"5. Converting bits to bytes:")
print(f"   Total bytes = {total_bits} bits / 8 bits/byte = {int(total_bytes)} bytes.")
print("-" * 20)

print(f"The smallest size of memory required is {int(total_bytes)} bytes.")
print(f"Final Equation: ceil( (24*60 + 39 + 35.244/60) ) * ceil( log2(1 + 2 + 2 + 4 + 2) ) / 8 = {int(total_bytes)}")

<<<740>>>