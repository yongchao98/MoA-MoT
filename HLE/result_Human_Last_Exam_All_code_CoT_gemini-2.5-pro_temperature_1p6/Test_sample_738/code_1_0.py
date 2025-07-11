import math

# Step 1: Calculate the total number of unique instructions.
# These are the different commands the rover can receive.
stop_cmd = 1
left_cmd = 2
right_cmd = 2
forward_cmd = 4
backward_cmd = 2
total_unique_instructions = stop_cmd + left_cmd + right_cmd + forward_cmd + backward_cmd

# Step 2: Calculate the minimum bits needed for one instruction.
# To represent N unique states, we need at least ceil(log2(N)) bits.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

# Step 3: Calculate the total number of instructions for one Martian day.
# A Martian sol is 24 hours, 39 minutes, and 35.244 seconds long.
# We need to find out how many Earth minutes this is.
sol_in_minutes = (24 * 60) + 39 + (35.244 / 60)

# Since an instruction is given every minute, we need to have enough instructions
# to cover the entire duration. We take the ceiling to ensure we cover the last partial minute.
total_instructions_needed = math.ceil(sol_in_minutes)

# Step 4: Calculate the total memory in bits.
total_bits = total_instructions_needed * bits_per_instruction

# Step 5: Convert the total bits to bytes (1 byte = 8 bits).
# We use math.ceil in case the total bits are not perfectly divisible by 8,
# though in this case they are.
total_bytes = math.ceil(total_bits / 8)

# Output the final equation with all the numbers.
print(f"The final calculation is:")
print(f"(ceil( (24*60 + 39 + 35.244/60) ) * ceil(log2({total_unique_instructions}))) / 8 = ?")
print(f"({total_instructions_needed} instructions * {bits_per_instruction} bits/instruction) / 8 bits/byte = {total_bytes} bytes")
print(f"\nThe smallest size of memory required is {total_bytes} bytes.")
