import re

def get_generators(group_string):
    """Calculates the number of generators for SU(N) or U(1) groups."""
    su_match = re.match(r"SU\((\d+)\)", group_string)
    u1_match = re.match(r"U\(1\)", group_string)
    
    if su_match:
        n = int(su_match.group(1))
        return n**2 - 1
    elif u1_match:
        return 1
    else:
        raise ValueError("Unsupported group format")

# --- Configuration ---
initial_group = "SU(3)"
final_groups = ["SU(2)", "U(1)"]

# --- Calculation ---
# 1. Calculate generators for the initial group G
num_gen_G = get_generators(initial_group)
n_g = int(re.search(r'\d+', initial_group).group())


# 2. Calculate generators for the final residual group H
num_gen_H = 0
final_group_calcs = []
for group in final_groups:
    num_gen_component = get_generators(group)
    num_gen_H += num_gen_component
    if "SU" in group:
      n_h = int(re.search(r'\d+', group).group())
      final_group_calcs.append(f"{n_h}^2 - 1")
    else:
      final_group_calcs.append("1")

# 3. Calculate the number of broken generators
broken_generators = num_gen_G - num_gen_H

# --- Output ---
print(f"Spontaneous symmetry breaking: {initial_group} -> {' x '.join(final_groups)}")
print("-" * 40)
print(f"Number of generators for initial group G ({initial_group}):")
print(f"dim(G) = {n_g}^2 - 1 = {num_gen_G}")
print("\nNumber of generators for residual group H:")
final_gen_calc_str = f"({') + ('.join(final_group_calcs)})"
print(f"dim(H) = {final_gen_calc_str} = {num_gen_H}")
print("\nNumber of broken generators is dim(G) - dim(H):")
print(f"{num_gen_G} - {num_gen_H} = {broken_generators}")

print("\nIn a non-Abelian gauge theory, the Higgs mechanism dictates that the number of")
print("massive gauge bosons is equal to the number of broken generators.")
print(f"Therefore, the unique condition resulting from this breaking is the existence of {broken_generators} massive gauge bosons.")