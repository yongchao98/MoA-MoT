import math

# This script calculates the smallest memory size in bytes to store one
# Martian day's worth of instructions for the Elon rover.

# Step 1: Define the number of states for each instruction type.
stop_instructions = 1
turn_left_instructions = 2
turn_right_instructions = 2
forward_instructions = 4
backward_instructions = 2

# Sum them to get the total number of distinct instructions.
total_distinct_instructions = (
    stop_instructions
    + turn_left_instructions
    + turn_right_instructions
    + forward_instructions
    + backward_instructions
)

# Step 2: Calculate the bits needed to store a single instruction.
# This is the ceiling of the base-2 logarithm of the total distinct instructions.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

# Step 3: Calculate the total number of instructions needed.
# A Martian day (sol) is 24 hours, 39 minutes, and 35.244 seconds.
# An instruction is issued each Earth minute.
mars_day_in_earth_minutes = (24 * 60) + 39 + (35.244 / 60)
# We need to cover the entire duration, so we take the ceiling.
total_instructions = math.ceil(mars_day_in_earth_minutes)

# Step 4: Calculate the total memory size.
# Total bits = total instructions * bits per instruction
# Total bytes = total bits / 8
total_bits = total_instructions * bits_per_instruction
total_bytes = total_bits / 8

print("--- Memory Calculation ---")
print(f"1. Total unique instructions: {stop_instructions} + {turn_left_instructions} + {turn_right_instructions} + {forward_instructions} + {backward_instructions} = {total_distinct_instructions}")
print(f"2. Minimum bits per instruction: ceil(log2({total_distinct_instructions})) = {bits_per_instruction}")
print(f"3. Total instructions for a Martian day (duration in minutes): ceil(24*60 + 39 + 35.244/60) = {total_instructions}")
print("\n--- Final Equation ---")
print(f"Memory in Bytes = ({total_instructions} instructions * {bits_per_instruction} bits/instruction) / 8 bits/byte")
print(f"Result = {int(total_bytes)} bytes")

<<<740>>>