# Step 1: Define the known variables based on the riddle's logic.
# From the simultaneous equations derived in the plan:
# 5 * L = 1400, where L is the number of lost pipes.
lost_pipes_numerator = 1400
lost_pipes_denominator = 5

# Step 2: Calculate the total number of lost pipes.
lost_pipes = lost_pipes_numerator / lost_pipes_denominator

# Step 3: The question asks how many pipes the tuner must find, which is half of the lost pipes.
divisor = 2
pipes_to_find = lost_pipes / divisor

# Step 4: Print the final equation with the calculated numbers.
print("The final calculation to find the number of pipes the tuner must realign is:")
print(f"{int(lost_pipes)} / {int(divisor)} = {int(pipes_to_find)}")