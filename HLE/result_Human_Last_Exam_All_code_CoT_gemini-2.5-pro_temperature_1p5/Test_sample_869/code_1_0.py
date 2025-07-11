import math

# Step-by-step calculation of the probability P = F/S

# Define constants
N_ITEMS = 25
N_INDIVIDUALS = 5
N_TYPES = 5
ITEMS_PER_TYPE = 5
ITEMS_PER_INDIVIDUAL = 5

# --- Step 1: Calculate S (Total number of distributions) ---
# The total number of ways to distribute 25 items to 5 individuals,
# with each receiving 5, is given by the multinomial coefficient:
# S = 25! / (5! * 5! * 5! * 5! * 5!) = 25! / (5!)^5
factorial_5 = math.factorial(ITEMS_PER_TYPE)
factorial_25 = math.factorial(N_ITEMS)
S = factorial_25 / (factorial_5**N_INDIVIDUALS)

# --- Step 2: Calculate F (Favorable number of distributions) ---
# A favorable distribution requires a unique dominant type for each individual.
# There are 5! ways to assign these dominant types to the individuals.
# We calculate the number of ways for one fixed assignment and multiply by 5!.

# For a fixed assignment, we count the ways based on the configuration matrix C.
# The number of ways to realize a configuration C is W(C) = (5!)^5 / product(c_ij!)

# Case A: Each individual gets 5 items of their dominant type.
# The matrix C is 5 times the identity matrix.
# c_ii=5, c_ij=0 for i!=j.
# product(c_ij!) = (5!)^5 * (0!)^20 = (5!)^5.
# W(C_A) = (5!)^5 / (5!)^5 = 1.
W_A = 1
n_A = 1 # There is only 1 such matrix for a fixed assignment.

# Case B: A mixed distribution where each individual has 2 of their dominant type.
# The matrix C has c_ii=2, and other entries are 0s and 1s.
# The number of such matrices equals the number of derangements of 5 items, D_5.
D_5 = 44
n_B = D_5
# For any such matrix, each column is a permutation of (2,1,1,1,0).
# The product of factorials of entries in one column is 2!*1!*1!*1!*0! = 2.
# For the whole matrix, product(c_ij!) = 2^5 = 32.
# W(C_B) = (5!)^5 / 32.
W_B = (factorial_5**N_TYPES) / (2**N_TYPES)

# Total favorable ways for a fixed assignment is the sum of ways for all valid matrices.
F_fixed_assignment = n_A * W_A + n_B * W_B

# Total F is F_fixed_assignment multiplied by the number of ways to assign dominant types (5!).
F = factorial_5 * F_fixed_assignment

# --- Step 3: Calculate the probability P = F/S ---
P = F / S

# --- Output the results ---
# The problem asks to output each number in the final equation.
# The final equation is P = F / S.
print("Step 1: Calculate the total number of distributions, S.")
print(f"S = {N_ITEMS}! / ({ITEMS_PER_TYPE}!) ^ {N_INDIVIDUALS}")
print(f"S = {S:.0f}")
print("-" * 20)

print("Step 2: Calculate the number of favorable distributions, F.")
print(f"F = {N_INDIVIDUALS}! * ( W_case_A + D_{N_INDIVIDUALS} * W_case_B )")
print(f"F = {factorial_5} * ( {n_A} * {W_A} + {n_B} * {W_B:.0f} )")
print(f"F = {F:.0f}")
print("-" * 20)

print("Step 3: Calculate the probability P = F / S.")
print(f"{P:.10f} = {F:.0f} / {S:.0f}")
