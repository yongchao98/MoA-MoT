# Define parameters based on the problem description
n_voters = 10
committee_size_k = 20
group_N_size = 8

# Step 1: Calculate s(N, W_1) for the core committee W_1
# From the analysis of the core property, to minimize satisfaction for group N,
# the committee W_1 is composed such that:
# - It contains the minimum required number of candidates from C_N = {1..8}, which is 7.
#   |W_1 intersect C_N| = 7
# - It contains the minimum required number of candidates from each P_i, which is 1.
#   |W_1 intersect P_i| = 1 for each i in N.
num_cands_from_C_N_for_W1 = 7
num_cands_from_each_P_i_for_W1 = 1
total_cands_from_all_P_i_for_W1 = group_N_size * num_cands_from_each_P_i_for_W1

# The satisfaction s(N, W) is defined as: 8 * |W intersect C_N| + sum(|W intersect P_i|)
s_N_W1 = group_N_size * num_cands_from_C_N_for_W1 + total_cands_from_all_P_i_for_W1

print("Calculation for W_1 (Core):")
print(f"s(N,W_1) = |N| * |W_1 intersect C_N| + sum(|W_1 intersect P_i| for i in N)")
print(f"s(N,W_1) = {group_N_size} * {num_cands_from_C_N_for_W1} + {total_cands_from_all_P_i_for_W1}")
print(f"s(N,W_1) = {s_N_W1}")
print("-" * 20)

# Step 2: Calculate s(N, W_2) for the EJR committee W_2
# From the analysis of EJR, to minimize satisfaction for group N,
# the committee W_2 is composed such that:
# - It contains 0 candidates from C_N to keep the main satisfaction term low.
#   |W_2 intersect C_N| = 0
# - To satisfy the EJR requirement that all voters in N are represented,
#   and to minimize cost, W_2 is constructed to have a total of 12 candidates
#   from the union of all P_i sets.
num_cands_from_C_N_for_W2 = 0
total_cands_from_all_P_i_for_W2 = 12

# The satisfaction s(N, W) is: 8 * |W intersect C_N| + sum(|W intersect P_i|)
s_N_W2 = group_N_size * num_cands_from_C_N_for_W2 + total_cands_from_all_P_i_for_W2

print("Calculation for W_2 (EJR):")
print(f"s(N,W_2) = |N| * |W_2 intersect C_N| + sum(|W_2 intersect P_i| for i in N)")
print(f"s(N,W_2) = {group_N_size} * {num_cands_from_C_N_for_W2} + {total_cands_from_all_P_i_for_W2}")
print(f"s(N,W_2) = {s_N_W2}")
print("-" * 20)

# Step 3: Compute the final ratio
if s_N_W2 != 0:
    ratio = s_N_W1 / s_N_W2
    print("Final Ratio Calculation:")
    print(f"Ratio = s(N,W_1) / s(N,W_2)")
    print(f"Ratio = {s_N_W1} / {s_N_W2}")
    print(f"Ratio = {ratio}")
else:
    print("s(N, W_2) is zero, ratio is undefined.")

# Final answer in specified format
# The fraction 64/12 simplifies to 16/3
final_answer = 16 / 3
print(f"\n<<<16/3>>>")
