import math

# Part 1: If n = 14, how many different exams can be created?
print("--- Part 1: Maximum exams for n=14 ---")
n1 = 14
k = 4
print(f"We are given n = {n1} total questions and k = {k} questions per exam.")
print("The condition is that any two exams can have at most one common question.")
print("\nStep 1: Find the maximum number of times a single question can appear.")
print("Let 'r' be the number of exams a specific question appears in.")
print("These 'r' exams all contain this one question. Each also contains k-1 = 3 other questions.")
print("Because any two of these exams can only share that one specific question, the sets of their other 3 questions must be disjoint.")
print("So, the total number of unique questions used by these 'r' exams is 1 (the common question) + r * (k-1).")
print("This total number of questions cannot exceed the total available questions, n.")
print("This gives us the inequality: 1 + r * (k-1) <= n")
print(f"Substituting the values: 1 + r * ({k}-1) <= {n1}")
print(f"1 + r * {k-1} <= {n1}")
print(f"r * {k-1} <= {n1 - 1}")
r_max = math.floor((n1 - 1) / (k - 1))
print(f"r <= {n1 - 1} / {k - 1} => r <= {round((n1 - 1) / (k - 1), 2)}")
print(f"Since r must be an integer, the maximum number of times any single question can appear is r_max = {r_max}.")

print("\nStep 2: Use this limit to find the maximum number of exams 'm'.")
print("The total number of question 'slots' across all 'm' exams is m * k.")
print("This total must equal the sum of the appearances of each question (let's call them r_i for each question i).")
print("So, m * k = sum(r_i for i=1 to n).")
print(f"Since each r_i must be less than or equal to r_max, the sum is at most n * r_max.")
print("This gives us a new inequality: m * k <= n * r_max")
print(f"Substituting the values: m * {k} <= {n1} * {r_max}")
m_max = math.floor((n1 * r_max) / k)
print(f"m * {k} <= {n1 * r_max}")
print(f"m <= {n1 * r_max} / {k}")
print(f"m <= {m_max}")
print(f"\nConclusion for Part 1: The maximum number of different exams that can be created is {m_max}.")

# Part 2: What is the minimum value of n needed to prepare 10 exams?
print("\n\n--- Part 2: Minimum n for 10 exams ---")
m2 = 10
k = 4
print(f"We are given m = {m2} exams and k = {k} questions per exam.")
print("We need to find the minimum possible value for n, the total number of questions.")

print("\nStep 1: Find a lower bound for n by counting pairs of questions.")
print("The number of distinct pairs of questions in a single exam is C(k, 2).")
pairs_per_exam = math.comb(k, 2)
print(f"C({k}, 2) = ({k} * {k-1}) / 2 = {pairs_per_exam}")
print(f"Across all {m2} exams, the total number of question pairs is {m2} * {pairs_per_exam} = {m2 * pairs_per_exam}.")
print("The condition 'at most one common question' implies that a specific pair of questions (e.g., Q1 and Q2) can appear together in at most one exam.")
print("The total number of possible pairs of questions from a pool of 'n' questions is C(n, 2).")
print("Therefore, the total pairs used in exams must be less than or equal to the total pairs available.")
print("This gives the inequality: m * C(k, 2) <= C(n, 2)")
inequality_rhs = m2 * pairs_per_exam * 2
print(f"{m2} * {pairs_per_exam} <= n * (n-1) / 2")
print(f"{inequality_rhs} <= n * (n-1)")

print("\nStep 2: Solve the inequality to find the minimum possible integer n.")
n_test = 1
while n_test * (n_test - 1) < inequality_rhs:
    n_test += 1
n_min_bound = n_test
print(f"We need to find the smallest integer n such that n*(n-1) >= {inequality_rhs}.")
print(f"For n = {n_min_bound - 1}, we have {n_min_bound - 1} * {n_min_bound - 2} = {(n_min_bound - 1) * (n_min_bound - 2)}, which is too small.")
print(f"For n = {n_min_bound}, we have {n_min_bound} * {n_min_bound - 1} = {n_min_bound * (n_min_bound - 1)}, which is large enough.")
print(f"This means n must be at least {n_min_bound}.")

print(f"\nStep 3: Check if n = {n_min_bound} is actually feasible.")
n_check = n_min_bound
print(f"Let's assume n = {n_check} and see if it leads to a contradiction.")
print("As in Part 1, let r_i be the number of exams question i appears in. We have the relation r_i <= (n-1)/(k-1).")
r_max_check = math.floor((n_check - 1) / (k - 1))
print(f"For n={n_check}, r_i <= ({n_check}-1)/({k}-1) = {round((n_check-1)/(k-1),2)}, so each question can appear in at most {r_max_check} exams.")
print("Now, let's count the total question slots: m * k = sum(r_i).")
total_slots = m2 * k
print(f"For m={m2}, the total slots are {m2} * {k} = {total_slots}.")
print(f"So, the sum of appearances of all {n_check} questions must be {total_slots}.")
print(f"However, the maximum possible sum of appearances is n * r_max = {n_check} * {r_max_check} = {n_check * r_max_check}.")
print(f"This leads to a contradiction, as {total_slots} cannot be greater than {n_check * r_max_check}.")
print(f"Therefore, n = {n_check} is not possible.")

print(f"\nStep 4: Conclude the minimum value of n.")
n_final = n_check + 1
print(f"Since n={n_check} is not possible, the minimum possible integer value for n is {n_final}.")
print(f"For n={n_final}, it is known that a set of exams can be constructed (specifically, a subset of the S(2, 4, 13) Steiner system).")
print(f"Conclusion for Part 2: The minimum value of n needed is {n_final}.")
