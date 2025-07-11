import math

# Step 1: Calculate the total number of unique instructions.
stop_instr = 1
left_instr = 2
right_instr = 2
forward_instr = 4
backward_instr = 2

total_unique_instructions = stop_instr + left_instr + right_instr + forward_instr + backward_instr
print(f"1. Total number of unique instructions:")
print(f"   Stop ({stop_instr}) + Turn Left ({left_instr}) + Turn Right ({right_instr}) + Forward ({forward_instr}) + Backward ({backward_instr}) = {total_unique_instructions}")

# Step 2: Calculate the minimum number of bits per instruction.
# We need ceil(log2(N)) bits to store N unique states.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))
print(f"\n2. Minimum bits per instruction:")
print(f"   ceil(log2({total_unique_instructions})) = {bits_per_instruction} bits")

# Step 3: Calculate the total number of instructions for one Martian day.
# A Martian day is 24 hours, 39 minutes, and 35 seconds.
# We convert this duration to minutes.
minutes_in_mars_day = (24 * 60) + 39 + (35 / 60)
# Since an instruction is sent EACH minute, we must cover the entire duration.
# So we need to take the ceiling to get the total number of instruction slots.
total_instructions_for_day = math.ceil(minutes_in_mars_day)
print(f"\n3. Total instructions for a Martian day:")
print(f"   Duration of a Martian day in minutes = {minutes_in_mars_day:.2f}")
print(f"   Total instructions needed = ceil({minutes_in_mars_day:.2f}) = {total_instructions_for_day}")

# Step 4: Calculate the total required memory in bytes.
total_bits = total_instructions_for_day * bits_per_instruction
# Convert bits to bytes (1 byte = 8 bits).
# Memory can't be a fraction of a byte, so we take the ceiling.
total_bytes = math.ceil(total_bits / 8)
print(f"\n4. Final memory calculation:")
print(f"   Total bits = {total_instructions_for_day} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
print(f"   Total bytes = {total_bits} bits / 8 bits/byte = {total_bytes} bytes")

print(f"\nThe smallest size of memory required is {total_bytes} bytes.")
print("Final equation:")
print(f"ceil( ( (24 * 60) + 39 + (35 / 60) ) ) * ceil( log2({stop_instr} + {left_instr} + {right_instr} + {forward_instr} + {backward_instr}) ) / 8 = {total_bytes} bytes")
print(f"Which simplifies to: {total_instructions_for_day} * {bits_per_instruction} / 8 = {total_bytes}")

