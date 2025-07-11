import math

# Step 1: Calculate the total number of distinct instructions.
stop_instr = 1
turn_left_instr = 2
turn_right_instr = 2
forward_instr = 4
backward_instr = 2

total_distinct_instructions = stop_instr + turn_left_instr + turn_right_instr + forward_instr + backward_instr

print("--- Step 1: Calculate the total number of distinct instructions ---")
print(f"The rover has the following instructions:")
print(f"Stop: {stop_instr} level")
print(f"Turn Left: {turn_left_instr} levels")
print(f"Turn Right: {turn_right_instr} levels")
print(f"Move Forward: {forward_instr} levels")
print(f"Move Backward: {backward_instr} levels")
print(f"Total distinct instructions = {stop_instr} + {turn_left_instr} + {turn_right_instr} + {forward_instr} + {backward_instr} = {total_distinct_instructions}")
print("-" * 20)

# Step 2: Determine the minimum number of bits required to represent these instructions.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

print("--- Step 2: Calculate the minimum bits per instruction ---")
print(f"To represent {total_distinct_instructions} unique states, we need to calculate ceil(log2({total_distinct_instructions})).")
print(f"The result is {bits_per_instruction} bits per instruction.")
print("-" * 20)

# Step 3: Calculate the total number of instructions for a whole day.
# Instructions are given each minute on Earth. An Earth day is used for calculation.
hours_in_day = 24
minutes_in_hour = 60
total_instructions_per_day = hours_in_day * minutes_in_hour

print("--- Step 3: Calculate total instructions for a day ---")
print(f"An instruction is sent every minute. An Earth day has {hours_in_day} hours and each hour has {minutes_in_hour} minutes.")
print(f"Total instructions = {hours_in_day} hours * {minutes_in_hour} minutes/hour = {total_instructions_per_day} instructions.")
print("-" * 20)

# Step 4: Calculate the total required memory in bytes.
total_bits = total_instructions_per_day * bits_per_instruction
# There are 8 bits in a byte.
bits_in_byte = 8
total_bytes = total_bits / bits_in_byte

print("--- Step 4: Calculate the total memory in bytes ---")
print(f"Total memory in bits = {total_instructions_per_day} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits.")
print(f"To convert bits to bytes, we divide by {bits_in_byte}.")
print(f"Total memory in bytes = {total_bits} bits / {bits_in_byte} bits/byte = {int(total_bytes)} bytes.")
print("-" * 20)

print(f"\nThe smallest size of memory required to store instructions for a whole day is {int(total_bytes)} bytes.")

print("<<<720>>>")