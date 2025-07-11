import math

# Step 1: Calculate the total number of distinct instruction types.
# Elon can stop (1), turn left (2 levels), turn right (2 levels),
# move forward (4 levels), or move backward (2 levels).
stop_instructions = 1
turn_left_instructions = 2
turn_right_instructions = 2
forward_instructions = 4
backward_instructions = 2

total_instruction_types = stop_instructions + turn_left_instructions + turn_right_instructions + forward_instructions + backward_instructions

# Step 2: Calculate the minimum bits needed to represent one instruction.
# This is the ceiling of log base 2 of the number of instruction types.
bits_per_instruction = math.ceil(math.log2(total_instruction_types))

# Step 3: Calculate the total number of instructions needed for one Martian day.
# A Martian day (sol) is 24 hours, 39 minutes, and 35.244 seconds long.
# We convert this duration to total Earth minutes.
# An instruction is issued each minute, so we need to cover the entire duration.
sol_hours = 24
sol_minutes = 39
sol_seconds = 35.244
total_minutes_in_sol = sol_hours * 60 + sol_minutes + sol_seconds / 60

# We need an instruction for each minute, including the last partial minute.
# So, we take the ceiling of the total minutes.
total_instructions_for_day = math.ceil(total_minutes_in_sol)

# Step 4: Calculate the total memory size in bytes.
# Total bits = bits per instruction * total instructions
# Total bytes = Total bits / 8 (rounded up)
total_bits_required = bits_per_instruction * total_instructions_for_day
total_bytes_required = math.ceil(total_bits_required / 8)

# Print the final equation with all the numbers.
print(f"The final calculation is: ceil( (ceil(log2({stop_instructions} + {turn_left_instructions} + {turn_right_instructions} + {forward_instructions} + {backward_instructions})) * ceil({sol_hours} * 60 + {sol_minutes} + {sol_seconds} / 60)) / 8 )")
print(f"= ceil( (ceil(log2({total_instruction_types})) * {total_instructions_for_day}) / 8 )")
print(f"= ceil( ({bits_per_instruction} * {total_instructions_for_day}) / 8 )")
print(f"= ceil( {total_bits_required} / 8 )")
print(f"= {int(total_bytes_required)} bytes")

# The final answer in the required format
# print(f"<<<{int(total_bytes_required)}>>>")