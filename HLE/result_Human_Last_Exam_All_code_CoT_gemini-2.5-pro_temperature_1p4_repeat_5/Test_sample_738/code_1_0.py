import math

# Step 1: Calculate the total number of unique instruction states.
# Stop: 1 state
# Turn Left: 2 speed levels
# Turn Right: 2 speed levels
# Move Forward: 4 speed levels
# Move Backward: 2 speed levels
stop_states = 1
turn_left_states = 2
turn_right_states = 2
forward_states = 4
backward_states = 2

total_states = stop_states + turn_left_states + turn_right_states + forward_states + backward_states

print(f"Calculating the total number of unique states:")
print(f"{stop_states} (stop) + {turn_left_states} (left) + {turn_right_states} (right) + {forward_states} (forward) + {backward_states} (backward) = {total_states} states")
print("-" * 30)

# Step 2: Calculate the minimum bits needed to represent one instruction.
# We need to find the smallest integer number of bits 'b' such that 2^b >= total_states.
# This is equivalent to ceil(log2(total_states)).
bits_per_instruction = math.ceil(math.log2(total_states))

print(f"Calculating the bits needed per instruction:")
print(f"ceil(log2({total_states})) = {bits_per_instruction} bits")
print("-" * 30)

# Step 3: Calculate the total number of instructions for one Martian day (Sol).
# A Martian Sol is 24 hours, 39 minutes, and 35.244 seconds.
# An instruction is sent each Earth minute.
sol_hours = 24
sol_minutes = 39
sol_seconds = 35.244
total_minutes_in_sol = (sol_hours * 60) + sol_minutes + (sol_seconds / 60)

# Since an instruction is sent each minute, we need to take the ceiling
# to ensure we have storage for the partial minute at the end.
total_instructions = math.ceil(total_minutes_in_sol)

print(f"Calculating the total instructions for a Martian day:")
print(f"Minutes in a Sol = ({sol_hours} * 60) + {sol_minutes} + ({sol_seconds} / 60) = {total_minutes_in_sol:.4f} minutes")
print(f"Total instructions needed = ceil({total_minutes_in_sol:.4f}) = {total_instructions}")
print("-" * 30)

# Step 4: Calculate total memory in bits.
total_bits = total_instructions * bits_per_instruction

print(f"Calculating the total bits required:")
print(f"{total_instructions} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
print("-" * 30)

# Step 5: Convert bits to bytes.
# There are 8 bits in a byte. We take the ceiling in case the total bits
# are not perfectly divisible by 8, as we cannot use a fraction of a byte.
bytes_in_total = math.ceil(total_bits / 8)

print(f"Converting total bits to bytes:")
print(f"{total_bits} bits / 8 bits/byte = {bytes_in_total} bytes")
print("-" * 30)
print(f"The smallest memory size required is {bytes_in_total} bytes.")
