import math

# Step 1: Calculate the total number of unique instructions.
stop_instructions = 1
turn_left_instructions = 2
turn_right_instructions = 2
forward_instructions = 4
backward_instructions = 2

total_distinct_instructions = (
    stop_instructions +
    turn_left_instructions +
    turn_right_instructions +
    forward_instructions +
    backward_instructions
)

# Step 2: Determine the number of bits required to represent all unique instructions.
# We need ceil(log2(N)) bits for N states.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

# Step 3: Calculate the total number of instructions for one Martian day.
# A Martian day is ~24.65 Earth hours. An instruction is given each Earth minute.
hours_in_martian_day = 24.65
minutes_per_hour = 60
total_instructions_in_a_day = int(hours_in_martian_day * minutes_per_hour)

# Step 4: Calculate the total memory required in bits and then convert to bytes.
# Total bits = total instructions * bits per instruction
total_bits = total_instructions_in_a_day * bits_per_instruction

# Total bytes = total bits / 8 bits per byte.
# We must use math.ceil to get the smallest whole number of bytes that can store the data.
total_bytes = math.ceil(total_bits / 8)

# Print the final equation with the numbers plugged in
print(f"To store enough instructions for a whole day on Mars, we need to solve the following:")
print(f"Total Unique Instructions = {stop_instructions} (stop) + {turn_left_instructions} (left) + {turn_right_instructions} (right) + {forward_instructions} (forward) + {backward_instructions} (backward) = {total_distinct_instructions}")
print(f"Bits per Instruction = ceil(log2({total_distinct_instructions})) = {bits_per_instruction} bits")
print(f"Total Instructions per Day = {hours_in_martian_day} hours * {minutes_per_hour} minutes/hour = {total_instructions_in_a_day} instructions")
print(f"Total Memory (bytes) = ceil( ({total_instructions_in_a_day} instructions * {bits_per_instruction} bits/instruction) / 8 bits/byte )")
print(f"Total Memory (bytes) = ceil( {total_bits} / 8 ) = {total_bytes} bytes")

print("\nThe smallest size of memory required is:")
print(total_bytes)
print("<<<740>>>")