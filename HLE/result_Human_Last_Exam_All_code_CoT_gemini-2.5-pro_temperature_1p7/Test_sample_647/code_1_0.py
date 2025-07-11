import math

# --- Problem Definition ---
# Number of voters
num_voters = 10
# Committee size
committee_size = 20
# Group N of voters
group_N = set(range(1, 9))

# --- Analysis of W1 (Core) ---
print("Step 1: Calculating the satisfaction s(N, W1) for the core committee.")
# Based on stability and core properties, W1 must contain the 8 common candidates.
num_common_candidates = 8
# To minimize s(N, W1), W1 must be filled with as many non-N candidates as possible.
# Candidates for voter 9 and 10 are used.
num_candidates_voter_9 = 4
num_candidates_voter_10 = 4

# The number of candidates in W1 from C_unique_N is the remainder.
num_unique_N_in_W1 = committee_size - num_common_candidates - num_candidates_voter_9 - num_candidates_voter_10

# s(N, W1) = |N| * (candidates from C_common) + (candidates from C_unique_N)
s_N_W1 = len(group_N) * num_common_candidates + num_unique_N_in_W1

print(f"A core committee W1 that minimizes satisfaction for N consists of:")
print(f"- {num_common_candidates} common candidates for group N.")
print(f"- {num_candidates_voter_9} candidates for voter 9.")
print(f"- {num_candidates_voter_10} candidates for voter 10.")
print(f"- {num_unique_N_in_W1} candidates from the unique pool for group N.")
print(f"The satisfaction s(N, W1) is calculated as |N|*|C_common intersect W1| + |C_unique_N intersect W1|.")
print(f"s(N, W1) = {len(group_N)} * {num_common_candidates} + {num_unique_N_in_W1} = {s_N_W1}")
print("-" * 20)

# --- Analysis of W2 (EJR) ---
print("Step 2: Calculating the satisfaction s(N, W2) for the EJR committee.")
# EJR implies W2 must contain all 8 common candidates for group N.
# EJR implies W2 must contain at least 2 candidates for voter 9 and 2 for voter 10.
# To minimize s(N, W2), we must minimize candidates from C_unique_N.
# This is achieved by maximizing candidates from voters 9 and 10.
max_cand_voter_9 = 4
max_cand_voter_10 = 4

# The number of candidates in W2 from C_unique_N is the remainder.
num_unique_N_in_W2 = committee_size - num_common_candidates - max_cand_voter_9 - max_cand_voter_10

# s(N, W2) = |N| * (candidates from C_common) + (candidates from C_unique_N)
s_N_W2 = len(group_N) * num_common_candidates + num_unique_N_in_W2

print(f"An EJR committee W2 that minimizes satisfaction for N must contain:")
print(f"- {num_common_candidates} common candidates due to the EJR rule for group N.")
print(f"- To minimize N's satisfaction, it includes the max possible from other voters: {max_cand_voter_9} for voter 9 and {max_cand_voter_10} for voter 10.")
print(f"- The remaining seats are filled by {num_unique_N_in_W2} candidates from N's unique pool.")
print(f"The satisfaction s(N, W2) is calculated as |N|*|C_common intersect W2| + |C_unique_N intersect W2|.")
print(f"s(N, W2) = {len(group_N)} * {num_common_candidates} + {num_unique_N_in_W2} = {s_N_W2}")
print("-" * 20)

# --- Final Calculation ---
print("Step 3: Calculating the final ratio.")
if s_N_W2 != 0:
    ratio = s_N_W1 / s_N_W2
    print(f"The final equation is s(N,W1) / s(N,W2).")
    print(f"The ratio is {s_N_W1} / {s_N_W2} = {ratio}")
else:
    print("Cannot calculate ratio because s(N, W2) is zero.")
