import math

# Step 1: Explain the derivation of the relationship between a1 and d.
print("Step 1: Determine the relationship between the first term (a1) and common difference (d) of sequence {a_n}.")
print("The sequence {a_n} is an arithmetic sequence, so a_n = a1 + (n-1)*d.")
print("The sequence {b_n} is defined as b_n = n*(n+1) / a_n.")
print("Since {b_n} is also an arithmetic sequence, the difference between consecutive terms is constant.")
print("This implies b_2 - b_1 = b_3 - b_2.")
print("Substituting the terms: b_1 = 2/a1, b_2 = 6/(a1+d), b_3 = 12/(a1+2d).")
print("The equation 6/(a1+d) - 2/a1 = 12/(a1+2d) - 6/(a1+d) simplifies to:")
print("a1^2 - 3*a1*d + 2*d^2 = 0")
print("Factoring this equation gives (a1 - d)*(a1 - 2d) = 0.")
print("This leads to two possible cases: Case 1 (a1 = d) and Case 2 (a1 = 2d).")
print("-" * 30)

# Step 2: Analyze Case 1 (a1 = d).
print("\nStep 2: Analyze Case 1: a1 = d")
print("If a1 = d, then a_n = d + (n-1)d = n*d.")
print("And b_n = n*(n+1) / (n*d) = (n+1)/d.")
print("Now, we find the sums S_99 and T_99.")
print("S_99 = sum_{k=1 to 99} (k*d) = d * (99 * 100 / 2) = 4950*d.")
print("T_99 = sum_{k=1 to 99} ((k+1)/d) = (1/d) * (sum_{k=1 to 99} k + sum_{k=1 to 99} 1) = (1/d) * (4950 + 99) = 5049/d.")
print("\nUsing the condition S_99 - T_99 = 99:")
print("4950*d - 5049/d = 99")
print("This simplifies to the quadratic equation: 50*d^2 - d - 51 = 0.")
# Solve 50*d^2 - d - 51 = 0
# Using quadratic formula: d = (-b +/- sqrt(b^2 - 4ac)) / 2a
a, b, c = 50, -1, -51
delta = b**2 - 4*a*c
d1 = (-b + math.sqrt(delta)) / (2*a)
d2 = (-b - math.sqrt(delta)) / (2*a)
print(f"Solving the equation gives d = {d1} or d = {d2}.")
print("Since the problem states d > 1, the valid solution for this case is d = 1.02.")
print("-" * 30)

# Step 3: Analyze Case 2 (a1 = 2d).
print("\nStep 3: Analyze Case 2: a1 = 2d")
print("If a1 = 2d, then a_n = 2d + (n-1)d = (n+1)*d.")
print("And b_n = n*(n+1) / ((n+1)*d) = n/d.")
print("Now, we find the sums S_99 and T_99.")
print("S_99 = sum_{k=1 to 99} ((k+1)*d) = d * (4950 + 99) = 5049*d.")
print("T_99 = sum_{k=1 to 99} (k/d) = (1/d) * (99 * 100 / 2) = 4950/d.")
print("\nUsing the condition S_99 - T_99 = 99:")
print("5049*d - 4950/d = 99")
print("This simplifies to the quadratic equation: 51*d^2 - d - 50 = 0.")
# Solve 51*d^2 - d - 50 = 0
a, b, c = 51, -1, -50
delta = b**2 - 4*a*c
d3 = (-b + math.sqrt(delta)) / (2*a)
d4 = (-b - math.sqrt(delta)) / (2*a)
print(f"Solving the equation gives d = {d3} or d = {round(d4, 4)}.")
print("Neither of these solutions satisfy the condition d > 1.")
print("-" * 30)

# Step 4: Final Conclusion and Verification
print("\nStep 4: Conclusion")
print("The only case that yields a valid result is Case 1, where d = 1.02.")
final_d = 1.02
s99_final = 4950 * final_d
t99_final = 5049 / final_d
print("\nFinal verification with d = 1.02:")
print(f"S_99 = 4950 * {final_d} = {s99_final}")
print(f"T_99 = 5049 / {final_d} = {t99_final}")
print(f"The final equation is S_99 - T_99 = {s99_final} - {t99_final} = {s99_final - t99_final}")
print("\nThis confirms the condition S_99 - T_99 = 99 is met.")
print(f"The value of d is {final_d}.")
