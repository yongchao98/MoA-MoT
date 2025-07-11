import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

# Step 1: Determine the constraints on the total number of pipes (T).

# Constraint from the fractions of the total pipes mentioned (1/3 and 2/5).
# T must be a multiple of lcm(3, 5).
constraint1 = lcm(3, 5)

# Constraint from the breakdown of the "lost pitch" pipes.
# This group's size is T/3.
# The group is broken into 1/4ths and 3/7ths, so its size must be a multiple of lcm(4, 7).
sub_constraint = lcm(4, 7)
# Therefore, T/3 must be a multiple of 28, so T must be a multiple of 3 * 28.
constraint2 = 3 * sub_constraint

# To satisfy all constraints, T must be a multiple of the lcm of all constraints.
total_pipes = lcm(constraint1, constraint2)

# Step 2: Use the known number of in-tune pipes to find the number of out-of-tune pipes.
in_tune_pipes = 200
total_lost_pipes = total_pipes - in_tune_pipes

# Step 3: Verify this against the poem's description of the out-of-tune groups.
# This step confirms our total_pipes calculation is consistent with all facts.
group_A_lost = total_pipes // 3 # "One-third of pipes fell out of tune"
group_B_affected = (2 * total_pipes) // 5 # "two-fifths caught the rising moon"
# The total lost pipes is the union of these groups. The overlap must be an integer.
# Overlap = |A| + |B| - |A U B|
overlap = group_A_lost + group_B_affected - total_lost_pipes
# We can see that with T=420, A=140, B=168, Lost=220, the overlap is 140+168-220 = 88, a valid number.

print(f"Based on the riddle's constraints, the total number of pipes in the cathedral is {total_pipes}.")
print(f"With {in_tune_pipes} pipes still singing pure, the number of 'lost' pipes is {total_pipes} - {in_tune_pipes} = {total_lost_pipes}.")

# Step 4: Calculate the final answer based on the last two lines.
# "How many must the tuner find / When just half the lost realign?"
# This means we need to find the number of pipes that do NOT realign.
divisor = 2
pipes_to_find = total_lost_pipes // divisor

print(f"\nThe poem states that half of the {total_lost_pipes} lost pipes realign.")
print(f"The number of pipes the tuner must find is the other half.")
print(f"\nThe final calculation is:")
print(f"{total_lost_pipes} / {divisor} = {pipes_to_find}")

<<<110>>>