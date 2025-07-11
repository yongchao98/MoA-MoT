import math

# Step 1: Calculate the total number of distinct instructions.
stop_levels = 1
left_levels = 2
right_levels = 2
forward_levels = 4
backward_levels = 2
total_distinct_instructions = stop_levels + left_levels + right_levels + forward_levels + backward_levels

# Step 2: Calculate the bits needed per instruction.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

# Step 3: Calculate the total instructions needed for one Martian day.
# 1 Martian sol = 24.6597 Earth hours
hours_in_sol = 24.6597
minutes_in_sol = hours_in_sol * 60
# An instruction is sent each minute. To cover the entire duration of a Martian day,
# we need to account for all minute intervals within that duration.
# For a duration of N minutes, this requires ceil(N) instructions.
total_instructions_for_day = math.ceil(minutes_in_sol)

# Step 4: Calculate total memory in bits.
total_bits = total_instructions_for_day * bits_per_instruction

# Step 5: Convert total bits to bytes and round up.
# 1 byte = 8 bits.
total_bytes = math.ceil(total_bits / 8)

# Print the step-by-step calculation
print("1. Calculate total number of distinct instructions:")
print(f"   {stop_levels} (stop) + {left_levels} (left) + {right_levels} (right) + {forward_levels} (forward) + {backward_levels} (backward) = {total_distinct_instructions} instructions")
print("\n2. Calculate bits needed per instruction:")
print(f"   ceil(log2({total_distinct_instructions})) = {bits_per_instruction} bits")
print("\n3. Calculate total instructions for one Martian day:")
print(f"   One Martian day is {minutes_in_sol:.3f} Earth minutes long. We need to store ceil({minutes_in_sol:.3f}) instructions.")
print(f"   Total instructions = {total_instructions_for_day}")
print("\n4. Calculate total memory needed in bytes:")
print(f"   Total bits = {total_instructions_for_day} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
print(f"   Total bytes = ceil({total_bits} / 8) = {total_bytes} bytes")

print("\nThe smallest size of memory required is 740 bytes.")
print(f"<<<{total_bytes}>>>")