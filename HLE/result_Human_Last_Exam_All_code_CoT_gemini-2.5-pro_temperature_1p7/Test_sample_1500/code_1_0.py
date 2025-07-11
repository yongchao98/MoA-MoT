# Plan:
# 1. Set n for the group G = SU(n). We will use n=4 as an example.
# 2. State the claim from the question for this n. The claimed value for b_2 is n-1.
# 3. Introduce the counterexample: the coadjoint orbit corresponding to the complex projective space CP^{n-1}.
# 4. State the known second Betti number for this space, which is 1.
# 5. Compare the claimed value with the actual value for the counterexample and conclude.

n = 4

# The claim is that for any orbit in SU(n), the second Betti number is n-1.
b2_claimed = n - 1

# A well-known coadjoint orbit for SU(n) is the complex projective space CP^{n-1}.
# For our example n=4, this is CP^3.
# The second Betti number of CP^{m} is always 1 for m >= 1.
b2_actual_for_counterexample = 1

print(f"The user's question (b) proposes that for G = SU({n}), the second Betti number b_2 of any coadjoint orbit is always {n}-1.")
print("Let's check this claim.")
print(f"The proposed equation is: b_2 = n - 1")
print(f"For n = {n}, the claimed value is: {n} - 1 = {b2_claimed}")
print("\nHowever, consider the coadjoint orbit which is the complex projective space CP^{n-1} = CP^{3}.")
print(f"The second Betti number for this space, b_2(CP^{n-1}), is known to be {b2_actual_for_counterexample}.")
print(f"\nComparing the values:")
print(f"Claimed b_2 = {b2_claimed}")
print(f"Actual b_2 for counterexample = {b2_actual_for_counterexample}")

if b2_claimed == b2_actual_for_counterexample:
    print("\nThe claim holds for this case.")
else:
    print(f"\nSince {b2_claimed} != {b2_actual_for_counterexample}, we have found a counterexample.")
    print("Therefore, the statement that the second Betti number is *always* n-1 is false.")
