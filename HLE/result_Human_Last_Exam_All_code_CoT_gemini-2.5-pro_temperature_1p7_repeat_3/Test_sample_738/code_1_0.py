import math

# Step 1: Calculate the total number of instructions for one Mars day (sol).
# A sol is 24 hours, 39 minutes, and 35.244 seconds in Earth time.
# An instruction is given each Earth minute.
hours_in_sol = 24
minutes_in_sol = 39
seconds_in_sol = 35.244

# Convert the duration to total minutes.
total_minutes = (hours_in_sol * 60) + minutes_in_sol + (seconds_in_sol / 60)
# Round up because an instruction is needed for the final partial minute.
total_instructions = math.ceil(total_minutes)

# Step 2: Determine the total number of distinct instructions.
stop_cmd = 1
left_cmd = 2
right_cmd = 2
forward_cmd = 4
backward_cmd = 2
total_distinct_instructions = stop_cmd + left_cmd + right_cmd + forward_cmd + backward_cmd

# Step 3: Calculate the minimum bits needed to represent one instruction.
# To represent N states, we need ceil(log2(N)) bits.
bits_per_instruction = math.ceil(math.log2(total_distinct_instructions))

# Step 4: Calculate total memory in bits, then convert to bytes.
total_bits = total_instructions * bits_per_instruction
total_bytes = total_bits / 8

# --- Final Calculation Output ---
print("Calculation Steps:")
print(f"1. Total instructions for a Mars day = ceil(({hours_in_sol} * 60) + {minutes_in_sol} + ({seconds_in_sol} / 60)) = {total_instructions}")
print(f"2. Number of unique instruction types = {stop_cmd} + {left_cmd} + {right_cmd} + {forward_cmd} + {backward_cmd} = {total_distinct_instructions}")
print(f"3. Bits needed per instruction = ceil(log2({total_distinct_instructions})) = {bits_per_instruction}")
print(f"4. Total bits needed = {total_instructions} instructions * {bits_per_instruction} bits/instruction = {total_bits} bits")
print(f"5. Total bytes needed = {total_bits} bits / 8 bits/byte = {int(total_bytes)} bytes")
