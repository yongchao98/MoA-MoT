import math

# Step 1: Calculate the total number of unique instructions
# Stop (1) + Turn Left (2) + Turn Right (2) + Forward (4) + Backward (2)
num_unique_instructions = 1 + 2 + 2 + 4 + 2

# Step 2: Determine the bits needed per instruction
# We need to find the smallest integer number of bits 'b' such that 2^b >= num_unique_instructions.
# This is equivalent to ceil(log2(num_unique_instructions)).
bits_per_instruction = math.ceil(math.log2(num_unique_instructions))

# Step 3: Calculate the total number of instructions for a day on Mars
# A Martian day (Sol) is 24 hours and 39 minutes in Earth time.
# An instruction is sent each Earth minute.
minutes_in_martian_day = (24 * 60) + 39

# Step 4: Calculate the total memory in bytes
# Total bits = minutes_in_martian_day * bits_per_instruction
# Total bytes = Total bits / 8. We use math.ceil because memory is allocated in whole bytes.
total_bits = minutes_in_martian_day * bits_per_instruction
total_bytes = math.ceil(total_bits / 8)

# Final Output
print(f"The number of unique instructions is: {num_unique_instructions}")
print(f"Bits needed per instruction: ceil(log2({num_unique_instructions})) = {bits_per_instruction}")
print(f"Minutes in a Martian day: (24 * 60) + 39 = {minutes_in_martian_day}")
print(f"\nTo find the total memory in bytes, we use the formula:")
print(f"Total Bytes = ceil((Total Instructions * Bits Per Instruction) / 8)")
print(f"Total Bytes = ceil(({minutes_in_martian_day} * {bits_per_instruction}) / 8)")
print(f"Total Bytes = ceil({total_bits} / 8)")
print(f"Total Bytes = ceil({total_bits / 8})")
print(f"Smallest memory size required: {total_bytes} bytes")
