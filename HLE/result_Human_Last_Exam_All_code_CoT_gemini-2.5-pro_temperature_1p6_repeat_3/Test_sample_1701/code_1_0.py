import math

# Step 1: Define Bloom filter parameters
m = 4  # number of bits in B
k = 2  # number of hash functions
n = 2  # number of inputs
N = k * n  # total number of hashings

print("Step-by-step calculation of the probability P(B'[i] = 1)")
print("-" * 50)
print(f"1. Parameters: m = {m} bits, k = {k} hashes, n = {n} items. Total hashes N = k*n = {N}.")

print("\n2. The goal is to find P(B'[i] = 1), where B'[i] is a bit in the derived filter.")
print("   This is P(B[j] != B[k]) for a corresponding pair of bits from B.")
print("   Let's calculate P(B[1] != B[3]).")

print("\n3. This probability can be expressed as:")
print("   P(B[1] != B[3]) = P(B[1]=0, B[3]=1) + P(B[1]=1, B[3]=0)")
print("   By symmetry, this simplifies to: 2 * P(B[1]=0, B[3]=1)")

print("\n4. To calculate P(B[1]=0, B[3]=1), we consider the locations of the N hashes:")
print("   - For B[1] to be 0, all hashes must avoid position 1. Prob of one hash avoiding 1 is (m-1)/m.")
prob_b1_zero_num = (m - 1) ** N
prob_b1_zero_den = m ** N
print(f"   P(B[1]=0) = (({m}-1)/{m})^{N} = (3/4)^{N} = {prob_b1_zero_num}/{prob_b1_zero_den}")

print("\n   - For B[1]=0 and B[3]=0, all hashes must avoid positions 1 and 3.")
print("     Prob of one hash avoiding 1 and 3 is (m-2)/m.")
prob_b1_b3_zero_num = (m - 2) ** N
prob_b1_b3_zero_den = m ** N
print(f"   P(B[1]=0, B[3]=0) = (({m}-2)/{m})^{N} = (2/4)^{N} = {prob_b1_b3_zero_num}/{prob_b1_b3_zero_den}")

print("\n   - The event (B[1]=0, B[3]=1) means all hashes avoid position 1, but not all avoid position 3.")
print("   P(B[1]=0, B[3]=1) = P(B[1]=0) - P(B[1]=0, B[3]=0)")
prob_p1_num = prob_b1_zero_num - prob_b1_b3_zero_num
prob_p1_den = prob_b1_zero_den
print(f"   P(B[1]=0, B[3]=1) = {prob_b1_zero_num}/{prob_b1_zero_den} - {prob_b1_b3_zero_num}/{prob_b1_b3_zero_den} = {prob_p1_num}/{prob_p1_den}")

print("\n5. Finally, combine the results for the final equation and answer.")
final_prob_num = 2 * prob_p1_num
final_prob_den = prob_p1_den
gcd_final = math.gcd(final_prob_num, final_prob_den)
final_simplified_num = final_prob_num // gcd_final
final_simplified_den = final_prob_den // gcd_final

print("\n--- FINAL EQUATION AND RESULT ---")
print(f"P(B'[i]=1) = 2 * [ P(B[1]=0) - P(B[1]=0, B[3]=0) ]")
print(f"P(B'[i]=1) = 2 * [ (({m}-1)/{m})^{N} - (({m}-2)/{m})^{N} ]")
print(f"P(B'[i]=1) = 2 * [ ({prob_b1_zero_num}/{prob_b1_zero_den}) - ({prob_b1_b3_zero_num}/{prob_b1_b3_zero_den}) ]")
print(f"P(B'[i]=1) = 2 * ( {prob_p1_num}/{prob_p1_den} )")
print(f"P(B'[i]=1) = {final_prob_num}/{final_prob_den}")
print(f"The final probability is {final_simplified_num}/{final_simplified_den}, which is {final_simplified_num/final_simplified_den}.")
