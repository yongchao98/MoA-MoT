def get_num_generators(group_name, n):
    """Calculates the number of generators for SU(N) or SO(N) Lie groups."""
    if group_name.upper() == "SU":
        return n**2 - 1
    elif group_name.upper() == "SO":
        return n * (n - 1) // 2
    else:
        return 0

# Initial group G = SU(3)
n_G = 3
dim_G = get_num_generators("SU", n_G)
print(f"The number of generators for the initial group G = SU({n_G}) is {dim_G}.")

# Residual group H = SO(3), assuming this was the intended group for option B
n_H = 3
dim_H = get_num_generators("SO", n_H)
print(f"The number of generators for the residual group H = SO({n_H}) is {dim_H}.")

# Calculate the number of broken generators
num_broken_generators = dim_G - dim_H

print("\nThe number of broken generators is the dimension of the initial group minus the dimension of the residual group.")
print("Final Equation:")
# Output each number in the final equation as requested
print(f"{dim_G} - {dim_H} = {num_broken_generators}")

print(f"\nThis corresponds to {num_broken_generators} broken generators, which matches option B.")
