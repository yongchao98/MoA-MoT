import math

# Step 1: Calculate the total number of unique instructions.
# The rover can stop, turn left/right, or move forward/backward.
stop_instr = 1
turn_left_instr = 2  # 2 speed levels
turn_right_instr = 2 # 2 speed levels
forward_instr = 4    # 4 speed levels
backward_instr = 2   # 2 speed levels

total_unique_instructions = (
    stop_instr + turn_left_instr + turn_right_instr +
    forward_instr + backward_instr
)

# Step 2: Determine the minimum number of bits required per instruction.
# We need to find the smallest integer 'n' such that 2^n is greater than or
# equal to the total number of unique instructions.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

# Step 3: Calculate the total number of instructions needed for one Martian day.
# A Martian day is 24 hours and 39 minutes long. An instruction is given
# every Earth minute.
minutes_in_mars_day = (24 * 60) + 39

# Step 4: Calculate the total memory required in bits.
total_bits = minutes_in_mars_day * bits_per_instruction

# Step 5: Convert the total bits to bytes.
# There are 8 bits in a byte. Since we cannot have a fraction of a byte for
# storage, we must round up to the nearest whole number.
total_bytes = math.ceil(total_bits / 8)

# Print the step-by-step calculation
print("--- Step 1: Count Unique Instructions ---")
print(f"Total unique instructions = {stop_instr} (stop) + {turn_left_instr} (left) + {turn_right_instr} (right) + {forward_instr} (forward) + {backward_instr} (backward) = {total_unique_instructions}")
print("\n--- Step 2: Calculate Bits per Instruction ---")
print(f"Minimum bits to represent {total_unique_instructions} instructions = ceil(log2({total_unique_instructions})) = {bits_per_instruction} bits")
print("\n--- Step 3: Calculate Total Instructions for a Martian Day ---")
print(f"Total minutes in a Martian day = (24 hours * 60 mins/hour) + 39 mins = {minutes_in_mars_day} minutes")
print(f"Total instructions needed = {minutes_in_mars_day}")
print("\n--- Step 4: Calculate Total Memory in Bits ---")
print(f"Total bits required = {minutes_in_mars_day} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
print("\n--- Step 5: Convert Bits to Bytes ---")
print(f"Total bytes required = ceil({total_bits} bits / 8 bits/byte) = {total_bytes} bytes")

print(f"\nThe smallest size of memory required is {total_bytes} bytes.")
<<<740>>>