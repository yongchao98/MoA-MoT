# All costs are in terms of multiplications (M) in the underlying field.
# We assume the cost of a squaring (S) is equal to one multiplication (M).

# Step 1: Compute Q = A + (-B), where A and B are affine points.
# The result Q is stored in extended coordinates.
# This is an affine-affine addition with the result in extended coordinates.
# According to the EFD formula 'add-2007-bl', this costs 9M.
cost_step1 = 9
print(f"Step 1: Compute Q = A - B. Inputs A, B are affine, output Q is extended.")
print(f"Cost of Step 1: {cost_step1}M")

# Step 2: Compute R = 2Q, where Q is in extended coordinates.
# This is a point doubling in extended coordinates.
# According to the EFD formula 'dbl-2007-bl', this costs 4S + 5M.
# Since S=M, the cost is 4M + 5M = 9M.
cost_step2 = 9
print(f"\nStep 2: Compute R = 2Q. Input Q is extended, output R is extended.")
print(f"Cost of Step 2: {cost_step2}M")

# Step 3: Compute P = R + (-B), where R is extended and -B is affine.
# This is a mixed-coordinate addition.
# According to the EFD formula 'madd-2007-bl', this costs 8M.
cost_step3 = 8
print(f"\nStep 3: Compute P = R - B. Input R is extended, input -B is affine, output P is extended.")
print(f"Cost of Step 3: {cost_step3}M")

# Calculate the total cost
total_cost = cost_step1 + cost_step2 + cost_step3
print(f"\nTotal cost for computing 2A - 3B using the strategy 2(A-B)-B is:")
print(f"{cost_step1}M (for A-B) + {cost_step2}M (for doubling) + {cost_step3}M (for final addition) = {total_cost}M")
print("\nThe final numerical answer for the cost is:")
print(total_cost)
