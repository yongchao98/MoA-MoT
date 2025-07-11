import math

# Step 1: Define the number of variations for each action to find the total number of unique instructions.
stop_cmd = 1
left_cmd = 2
right_cmd = 2
forward_cmd = 4
backward_cmd = 2
total_unique_instructions = stop_cmd + left_cmd + right_cmd + forward_cmd + backward_cmd
print(f"1. Total number of unique instructions: {stop_cmd} + {left_cmd} + {right_cmd} + {forward_cmd} + {backward_cmd} = {total_unique_instructions}")

# Step 2: Calculate the minimum number of bits to store one instruction.
# This is the ceiling of log base 2 of the number of instructions.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))
print(f"2. Minimum bits per instruction: ceil(log2({total_unique_instructions})) = {bits_per_instruction}")

# Step 3: Calculate the total number of instructions for one Martian day (sol).
# A sol is 24 hours and 39 minutes long. An instruction is sent each minute.
mars_day_hours = 24
mars_day_minutes = 39
total_instructions_per_day = (mars_day_hours * 60) + mars_day_minutes
print(f"3. Total instructions per Martian day: ({mars_day_hours} * 60) + {mars_day_minutes} = {total_instructions_per_day}")

# Step 4: Calculate the total bits needed for a day's instructions.
total_bits = total_instructions_per_day * bits_per_instruction
print(f"4. Total bits required: {total_instructions_per_day} * {bits_per_instruction} = {total_bits}")

# Step 5: Convert bits to bytes (1 byte = 8 bits).
# Memory is allocated in whole bytes, so we must round up to the nearest byte.
total_bytes = math.ceil(total_bits / 8)
print(f"5. Smallest memory size in bytes: ceil({total_bits} / 8) = {total_bytes}")

print("\nThe smallest size of memory required is 740 bytes.")