import math

# Step 1: Calculate the total number of unique instructions.
stop_levels = 1
left_turn_levels = 2
right_turn_levels = 2
forward_levels = 4
backward_levels = 2

total_instruction_types = stop_levels + left_turn_levels + right_turn_levels + forward_levels + backward_levels

# Step 2: Calculate the minimum number of bits to store one instruction.
# We need to find the smallest integer 'b' such that 2^b >= total_instruction_types.
# This is equivalent to ceil(log2(total_instruction_types)).
bits_per_instruction = math.ceil(math.log2(total_instruction_types))

# Step 3: Calculate the total number of instructions needed for one Martian day (sol).
# A sol is 24h 39m 35.244s. An instruction is sent each Earth minute.
sol_hours = 24
sol_minutes = 39
sol_seconds = 35.244

minutes_in_sol = (sol_hours * 60) + sol_minutes + (sol_seconds / 60)
# We need to cover every minute, so we take the ceiling.
total_instructions_per_day = math.ceil(minutes_in_sol)

# Step 4: Calculate the total memory in bits.
total_bits = total_instructions_per_day * bits_per_instruction

# Step 5: Convert bits to bytes (1 byte = 8 bits).
total_bytes = total_bits / 8

# Print the final result and the equation used to find it.
# The number of bytes must be an integer.
total_bytes = int(total_bytes)

print(f"To store enough instructions for a whole day on Mars, we need:")
print(f"({total_instructions_per_day} instructions * {bits_per_instruction} bits/instruction) / 8 bits/byte = {total_bytes} bytes")