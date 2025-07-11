import math

# --- Constants from the problem statement ---
k = 100  # Committee size

# --- Part 1: Finding s1 for Proportional Justified Representation (PJR) ---
print("--- Finding s1 for Proportional Justified Representation (PJR) ---")
print(f"The committee size k is {k}.")
print("We want the smallest number of voters n such that there is a committee W satisfying PJR, while leaving voter 1 unsatisfied.")
print("Voter 1's ballot is A(1) = {a, b, c, x}. Being unsatisfied means W ∩ {a, b, c, x} = ∅.")
print("\nPJR states: For any l-cohesive group N' with |N'| >= n/k * l, at least one voter in N' must be satisfied.")
print("To find the smallest n, we create a 'best-case' scenario for PJR by choosing a W that satisfies as many other voters as possible.")
print("The ballots for voters 2-6 are A(2,3)={a,b,c,y} and A(4,5,6)={a,b,c,z}.")
print("By choosing W to contain {y, z}, we satisfy voters 2 through 6.")
print("We can then construct the rest of the profile and W such that only voter 1 is unsatisfied.")
print("\nIn this scenario, the only group of fully unsatisfied voters is N' = {1}.")
print("For PJR to hold, this group must NOT meet the size condition that would guarantee it representation.")
print("For N' = {1}, the ballot has 4 candidates, so the group is 4-cohesive (l=4). The group size is |N'|=1.")
l1 = 4
size_N1 = 1
print(f"The PJR size condition for a violation is |N'| >= n/k * l, which is: {size_N1} >= n / {k} * {l1}.")
print("To ensure PJR holds, this condition must be false:")
print(f"Equation: {size_N1} < n / {k} * {l1}")
print(f"Solving for n: {size_N1 * k} < {l1} * n")
print(f"Further solving: {size_N1 * k / l1} < n")
n_bound_s1 = size_N1 * k / l1
print(f"Resulting inequality: {n_bound_s1} < n")
s1 = math.floor(n_bound_s1) + 1
print(f"The smallest integer n satisfying this is {s1}.")
print(f"Therefore, s1 = {s1}")

# --- Part 2: Finding s2 for Extended Justified Representation (EJR) ---
print("\n\n--- Finding s2 for Extended Justified Representation (EJR) ---")
print("EJR states: For any l-cohesive group N' with |N'| >= n/k * l, |(∩ A(i)) ∩ W| >= l.")
print("As before, since voter 1 is unsatisfied, we have W ∩ {a, b, c, x} = ∅.")
print("\nConsider the group N' = {1, 2, 3, 4, 5, 6}. Its size is |N'| = 6.")
size_N2 = 6
print("The intersection of their ballots is ∩A(i) = {a, b, c}. This group is 3-cohesive, so l = 3.")
l2 = 3
print("\nThe EJR rule for this group is: IF |N'| >= n/k * l, THEN |{{a, b, c}} ∩ W| >= l.")
print(f"The conclusion |{{a, b, c}} ∩ W| >= {l2} means W must contain all of {a, b, c}.")
print("This contradicts the fact that W ∩ {a, b, c, x} = ∅.")
print("\nTherefore, for the situation to be possible, the premise of the EJR rule must be false.")
print("The premise is |N'| >= n/k * l. We need this to be false, so: |N'| < n/k * l.")
print(f"Equation: {size_N2} < n / {k} * {l2}")
print(f"Solving for n: {size_N2 * k} < {l2} * n")
print(f"Further solving: {size_N2 * k / l2} < n")
n_bound_s2 = size_N2 * k / l2
print(f"Resulting inequality: {n_bound_s2} < n")
s2 = math.floor(n_bound_s2) + 1
print(f"The smallest integer n satisfying this is {s2}.")
print(f"Therefore, s2 = {s2}")

# --- Final Result ---
print("\n--- Final Answer ---")
final_answer = (s1, s2)
print(f"<<<{final_answer}>>>")