import math

# Step 1: Calculate the total number of unique instructions.
# The rover has several actions with different speed levels.
stop_instructions = 1
turn_left_instructions = 2
turn_right_instructions = 2
forward_instructions = 4
backward_instructions = 2

# Sum them up to get the total number of unique commands.
total_unique_instructions = (stop_instructions +
                             turn_left_instructions +
                             turn_right_instructions +
                             forward_instructions +
                             backward_instructions)

# Step 2: Determine the minimum number of bits to represent one instruction.
# To represent N unique states, we need ceil(log2(N)) bits.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

# Step 3: Calculate the total number of instructions for a whole Martian day.
# A Martian day (sol) is 24 hours, 39 minutes, and 35.244 seconds.
# We convert this duration into Earth minutes.
hours_in_sol = 24
minutes_in_sol = 39
seconds_in_sol = 35.244
total_minutes_in_sol = (hours_in_sol * 60) + minutes_in_sol + (seconds_in_sol / 60)

# An instruction is sent each Earth minute. To cover the entire Martian day,
# we round up the number of minutes to the next whole number.
total_instructions_for_a_day = math.ceil(total_minutes_in_sol)

# Step 4: Calculate the total memory size in bits.
total_bits = total_instructions_for_a_day * bits_per_instruction

# Step 5: Convert the total bits to bytes.
# There are 8 bits in a byte. Memory is allocated in whole bytes,
# so we take the ceiling of the division.
total_bytes = math.ceil(total_bits / 8)

# Print the final calculation and result.
print(f"To store {total_unique_instructions} unique instructions, we need {bits_per_instruction} bits per instruction.")
print(f"A Martian day requires {total_instructions_for_a_day} instructions (one per Earth minute).")
print("The calculation for the total memory in bytes is:")
print(f"({total_instructions_for_a_day} instructions * {bits_per_instruction} bits/instruction) / 8 bits/byte = {total_bytes} bytes")
print(f"The smallest size of memory required is {total_bytes} bytes.")
