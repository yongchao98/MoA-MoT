import math
from fractions import Fraction

# Step 1: Define the aFC values and their ranks. Using Fractions for precision.
# aFC ranks: 1/3 (1), 1/2 (2), 3/2 (3), 2 (4), 3 (5)
afc_values = [Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(2, 1), Fraction(3, 1)]
ranks = {val: i+1 for i, val in enumerate(sorted(afc_values))}

print("Goal: Find the sequence of SNPs (by rank) that satisfies the expression conditions.")
print("The conditions are that for some SNP order (f1..f5), there exist k1 and k2 in {2,3,4} such that:")
print("E(k1) = 2.0 and E(k2) = 5.0")
print("Where E(k) = f_k * ( product(f_{i}, i>k) + product(f_{j}, j<k) )\n")

solution_found = False
for p in sorted(list(set(afc_values))): # Use a fixed sorted order for permutations
    perm = list(p)
    f = perm # f_1=f[0], f_2=f[1], ...

    # Calculate expressions for M/M SNP at positions k=2, 3, 4
    # E(k=2) = f_2 * (f_3*f_4*f_5 + f_1)
    val_k2 = f[1] * (f[2]*f[3]*f[4] + f[0])

    # E(k=3) = f_3 * (f_4*f_5 + f_1*f_2)
    val_k3 = f[2] * (f[3]*f[4] + f[0]*f[1])

    # E(k=4) = f_4 * (f_5 + f_1*f_2*f_3)
    val_k4 = f[3] * (f[4] + f[0]*f[1]*f[2])

    results = [val_k2, val_k3, val_k4]

    # Check if the results contain both 2 and 5
    found_2 = any(math.isclose(float(r), 2.0) for r in results)
    found_5 = any(math.isclose(float(r), 5.0) for r in results)

    if found_2 and found_5:
        solution_ranks = [ranks[val] for val in f]
        
        print("Solution Found!")
        print(f"The order of aFC values is: {', '.join(map(str, f))}")
        print(f"The corresponding rank order is: {''.join(map(str, solution_ranks))}\n")

        print("Verification:")
        print(f"Expression if M/M is at pos 2: {float(val_k2):.2f}")
        print(f"Expression if M/M is at pos 3: {float(val_k3):.2f}")
        print(f"Expression if M/M is at pos 4: {float(val_k4):.2f}")
        
        # Format the final answer string as requested
        final_ranks = solution_ranks
        if final_ranks[0] > final_ranks[-1]:
            final_ranks.reverse()
        final_answer_str = "".join(map(str, final_ranks))

        print(f"\nThe SNP order by rank is {final_answer_str}")
        solution_found = True
        break

if not solution_found:
    print("No solution found that perfectly matches the criteria {2.0, 5.0}.")

<<<34125>>>