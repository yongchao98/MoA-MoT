import math

# Step 1: Calculate the total number of unique instructions.
stop_levels = 1
turn_left_levels = 2
turn_right_levels = 2
forward_levels = 4
backward_levels = 2

total_unique_instructions = stop_levels + turn_left_levels + turn_right_levels + forward_levels + backward_levels
print(f"1. Total number of unique instructions:")
print(f"   Stop ({stop_levels}) + Turn Left ({turn_left_levels}) + Turn Right ({turn_right_levels}) + Forward ({forward_levels}) + Backward ({backward_levels}) = {total_unique_instructions}")
print("-" * 30)

# Step 2: Determine the bits needed per instruction.
# To represent N states, we need ceil(log2(N)) bits.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))
print(f"2. Minimum bits per instruction:")
print(f"   ceil(log2({total_unique_instructions})) = {bits_per_instruction} bits")
print("-" * 30)

# Step 3: Calculate the total number of instructions for one Martian day.
# A Martian day (sol) is approximately 24 hours, 39.5 minutes long in Earth time.
hours_in_sol = 24
minutes_in_sol_extra = 39.5
total_minutes_in_sol = (hours_in_sol * 60) + minutes_in_sol_extra
# An instruction is sent each minute, so we need to store instructions for the entire duration.
# We take the ceiling to account for the last partial minute.
total_instructions_for_day = math.ceil(total_minutes_in_sol)
print(f"3. Total instructions for a Martian day:")
print(f"   Duration of a Martian day in Earth minutes = ({hours_in_sol} * 60) + {minutes_in_sol_extra} = {total_minutes_in_sol}")
print(f"   Total instructions needed = ceil({total_minutes_in_sol}) = {total_instructions_for_day}")
print("-" * 30)

# Step 4: Calculate the total memory in bits.
total_bits = total_instructions_for_day * bits_per_instruction
print(f"4. Total memory required in bits:")
print(f"   {total_instructions_for_day} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
print("-" * 30)

# Step 5: Convert bits to bytes.
# 1 byte = 8 bits. The result must be an integer, as we are calculating the size of memory.
bits_in_a_byte = 8
total_bytes = total_bits / bits_in_a_byte
print(f"5. Total memory required in bytes:")
print(f"   {total_bits} bits / {bits_in_a_byte} bits/byte = {int(total_bytes)} bytes")
print("-" * 30)

# Final summary equation
print("Final Calculation:")
print(f"({bits_per_instruction} bits/instruction * {total_instructions_for_day} instructions) / {bits_in_a_byte} bits/byte = {int(total_bytes)} bytes")
<<<740>>>