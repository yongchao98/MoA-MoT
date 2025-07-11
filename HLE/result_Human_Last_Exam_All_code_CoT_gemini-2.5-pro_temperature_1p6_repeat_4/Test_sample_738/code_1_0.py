import math

# Step 1: Define the number of states for each action.
stop_instructions = 1
left_turn_instructions = 2
right_turn_instructions = 2
forward_instructions = 4
backward_instructions = 2

# Calculate the total number of distinct instructions.
total_distinct_instructions = (
    stop_instructions
    + left_turn_instructions
    + right_turn_instructions
    + forward_instructions
    + backward_instructions
)

# Step 2: Determine the minimum bits required per instruction.
# To represent N states, we need ceil(log2(N)) bits.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

# Step 3: Calculate the total number of instructions for one Martian day.
# A Martian day (sol) is 88775.244 seconds.
seconds_in_mars_day = 88775.244
# An instruction is sent each Earth minute (60 seconds).
seconds_per_instruction = 60
# Total instructions is the ceiling of the duration divided by the interval.
total_instructions_for_a_day = math.ceil(seconds_in_mars_day / seconds_per_instruction)

# Step 4: Calculate the total memory in bits and then in bytes.
total_bits = total_instructions_for_a_day * bits_per_instruction
total_bytes = total_bits / 8

# Print the breakdown of the calculation
print(f"1. Total unique instructions: {stop_instructions} + {left_turn_instructions} + {right_turn_instructions} + {forward_instructions} + {backward_instructions} = {total_distinct_instructions}")
print(f"2. Bits needed per instruction: ceil(log2({total_distinct_instructions})) = {bits_per_instruction}")
print(f"3. Total instructions for a Martian day: ceil({seconds_in_mars_day} / {seconds_per_instruction}) = {total_instructions_for_a_day}")
print(f"4. Total memory in bytes: ({bits_per_instruction} * {total_instructions_for_a_day}) / 8 = {int(total_bytes)}")
print("\nFinal equation and result:")
# The final equation with all numbers included
print(f"(ceil(log2({stop_instructions} + {left_turn_instructions} + {right_turn_instructions} + {forward_instructions} + {backward_instructions})) * ceil({seconds_in_mars_day} / {seconds_per_instruction})) / 8 = {int(total_bytes)}")