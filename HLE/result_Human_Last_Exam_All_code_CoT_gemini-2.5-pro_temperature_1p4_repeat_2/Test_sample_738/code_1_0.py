import math

# Step 1: Calculate the total number of unique instructions.
# The rover's CPU can issue instructions for stop, turn left, turn right, forward, or backward.
# Each movement command has a specified number of speed levels.
num_stop_instructions = 1
num_turn_left_instructions = 2
num_turn_right_instructions = 2
num_forward_instructions = 4
num_backward_instructions = 2

# Summing them up to get the total number of unique states.
total_unique_instructions = (
    num_stop_instructions +
    num_turn_left_instructions +
    num_turn_right_instructions +
    num_forward_instructions +
    num_backward_instructions
)

print(f"Step 1: Calculate the total number of unique instructions.")
print(f"Total unique instructions = {num_stop_instructions} (stop) + {num_turn_left_instructions} (turn left) + {num_turn_right_instructions} (turn right) + {num_forward_instructions} (forward) + {num_backward_instructions} (backward)")
print(f"Total unique instructions = {total_unique_instructions}\n")


# Step 2: Calculate the minimum bits required to store one instruction.
# To represent N unique states, we need to find the smallest integer 'b' where 2^b >= N.
# This is equivalent to the ceiling of log base 2 of N.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

print(f"Step 2: Calculate bits per instruction.")
print(f"To represent {total_unique_instructions} states, we need ceil(log2({total_unique_instructions})) bits.")
print(f"Bits per instruction = {bits_per_instruction}\n")


# Step 3: Calculate the total number of instructions for one Martian day (sol).
# A sol's duration is ~24 hours, 39 minutes, and 35 seconds in Earth time.
# An instruction is given each Earth minute.
hours_in_sol = 24
minutes_part_of_sol = 39
seconds_part_of_sol = 35

total_minutes_in_sol = (hours_in_sol * 60) + minutes_part_of_sol + (seconds_part_of_sol / 60)

# We must cover the entire Martian day, so we take the ceiling of the total minutes.
total_instructions_for_a_day = math.ceil(total_minutes_in_sol)

print(f"Step 3: Calculate total instructions for a Martian day.")
print(f"A Martian day in Earth time is {hours_in_sol}h {minutes_part_of_sol}m {seconds_part_of_sol}s, or {total_minutes_in_sol:.2f} minutes.")
print(f"Total instructions needed = ceil({total_minutes_in_sol:.2f}) = {total_instructions_for_a_day}\n")


# Step 4: Calculate the total memory size in bytes.
# Total memory (bits) = Total Instructions * Bits per Instruction
# Total memory (bytes) = Total memory (bits) / 8
bits_in_a_byte = 8
total_memory_bits = total_instructions_for_a_day * bits_per_instruction
total_memory_bytes = total_memory_bits / bits_in_a_byte

print(f"Step 4: Calculate the final memory size in bytes.")
print(f"Total bits = {total_instructions_for_a_day} instructions * {bits_per_instruction} bits/instruction = {total_memory_bits} bits.")
print(f"Total memory in bytes = {total_memory_bits} bits / {bits_in_a_byte} bits/byte")
print(f"The smallest size of memory required is: {int(total_memory_bytes)} bytes.")
