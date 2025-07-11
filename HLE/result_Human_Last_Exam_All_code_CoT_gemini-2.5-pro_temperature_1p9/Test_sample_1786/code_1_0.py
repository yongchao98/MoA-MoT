# The final formula combines the zero and non-zero cases using the 'plus' connective (external choice).
# Let's represent the formula components as strings.

r = 'r'
z = 'z'
nz = 'nz'

# Formula for the zero-branch: transition to state 'z'.
# This simply produces the new state literal S_z.
zero_branch = f"S_{z}"

# Formula for the non-zero branch: decrement counter 'r' and transition to state 'nz'.
# This consumes a C_r resource and produces the new state literal S_nz.
# The linear implication 'multimap' is represented by '-o'.
nonzero_branch = f"(C_{r} -o S_{nz})"

# The full formula F(r,z,nz) is the external choice between these two branches.
# External choice 'plus' is represented by '⊕'.
F_formula = f"{zero_branch} ⊕ {nonzero_branch}"

print("The formula F(r,z,nz) is:")
# The linear logic symbols can be represented with ASCII or Unicode.
# Using Unicode for clarity.
# ⊸ is \u22b8
# ⊕ is \u2295
pretty_formula = F_formula.replace('-o', '\u22b8').replace('⊕', '\u2295')
print(pretty_formula)

final_answer = zero_branch + " ⊕ " + nonzero_branch

# Printing the formula again, as requested by the output format.
# Let's use characters available on most keyboards for the final output string.
# S_z represents the state atom for z
# C_r represents the counter atom for r
# -o represents the linear implication (multimap)
# (+) represents the additive disjunction (plus)

# Each component needs to be part of the final equation string
F_final_str = f"S_z (+) (C_r -o S_nz)"
print("Here is the formula written with basic characters:")
print("F(r,z,nz) = S_z (+) (C_r -o S_nz)")

final_answer_repr = "S_z \u2295 (C_r \u22b8 S_{nz})"