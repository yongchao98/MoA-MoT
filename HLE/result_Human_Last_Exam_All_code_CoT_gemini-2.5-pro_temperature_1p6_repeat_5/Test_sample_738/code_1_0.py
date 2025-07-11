import math

# Step 1: Calculate the total number of unique instructions.
# There are 5 types of actions with different numbers of levels.
stop_instr = 1
turn_left_instr = 2
turn_right_instr = 2
forward_instr = 4
backward_instr = 2

# Summing them up to get the total number of possible unique instructions.
total_unique_instructions = stop_instr + turn_left_instr + turn_right_instr + forward_instr + backward_instr

print("--- Step 1: Calculate Unique Instructions ---")
print(f"The number of unique instructions is the sum of all action variations:")
print(f"Equation: {stop_instr} (stop) + {turn_left_instr} (left) + {turn_right_instr} (right) + {forward_instr} (forward) + {backward_instr} (backward) = {total_unique_instructions}\n")


# Step 2: Calculate the minimum bits needed to represent one instruction.
# To represent N states, we need ceil(log2(N)) bits.
bits_per_instruction = math.ceil(math.log2(total_unique_instructions))

print("--- Step 2: Calculate Bits per Instruction ---")
print(f"To represent {total_unique_instructions} unique instructions, we need a certain number of bits.")
print(f"Equation: ceil(log2({total_unique_instructions})) = {bits_per_instruction} bits\n")


# Step 3: Calculate the total number of instructions for one Martian day.
# A Martian day (sol) is 24 hours and 39 Earth minutes.
# An instruction is given each Earth minute.
hours_in_sol = 24
minutes_in_sol = 39
total_minutes_in_sol = (hours_in_sol * 60) + minutes_in_sol

print("--- Step 3: Calculate Total Instructions for a Martian Day ---")
print("A Martian day is 24 hours and 39 minutes. An instruction is sent every minute.")
print(f"Equation: ({hours_in_sol} * 60) + {minutes_in_sol} = {total_minutes_in_sol} total instructions\n")


# Step 4: Calculate the total memory required in bits.
total_bits = total_minutes_in_sol * bits_per_instruction

print("--- Step 4: Calculate Total Memory in Bits ---")
print(f"Total bits are the total instructions multiplied by the bits per instruction.")
print(f"Equation: {total_minutes_in_sol} * {bits_per_instruction} = {total_bits} bits\n")


# Step 5: Convert the total bits to bytes.
# There are 8 bits in a byte. We must round up for storage (ceiling).
bits_per_byte = 8
total_bytes = math.ceil(total_bits / bits_per_byte)

print("--- Step 5: Convert Bits to Bytes ---")
print(f"To get the size in bytes, we divide the total bits by {bits_per_byte} and round up.")
print(f"Equation: ceil({total_bits} / {bits_per_byte}) = {int(total_bytes)}")
print("\nFinal Answer:")
print(f"The smallest size of memory required is {int(total_bytes)} bytes.")
