# This script calculates the infimum of c for which a given Markov chain is transient.

# Plan:
# 1. For a random walk on the integers, its long-term behavior (transience or recurrence)
#    is determined by its asymptotic drift and variance.
# 2. We calculate the drift (mean jump) mu_k = E[X_{n+1} - k | X_n = k].
# 3. We calculate the second moment of the jump sigma_k^2 = E[(X_{n+1} - k)^2 | X_n = k].
# 4. We apply Lamperti's criterion, which states that the chain is transient if the exponent p
#    in an associated p-series is greater than 1. This exponent is determined by the
#    asymptotic behavior of (2 * mu_k * k) / sigma_k^2.
# 5. We solve the inequality p > 1 for c to find the condition for transience.
# 6. The infimum is the boundary value of this condition.

print("Step 1: Calculate the asymptotic drift mu_k")
# The drift is the expected change in position from state k.
# mu_k = (-2)*P_{k,k-2} + (-1)*P_{k,k-1} + (1)*P_{k,k+1} + (2)*P_{k,k+2}
# mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)
# mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2 = 2*c/k
drift_coeff_A = 2
print(f"For large k, the drift mu_k is of the form A/k, where A = {drift_coeff_A}*c.")
print("-" * 40)

print("Step 2: Calculate the asymptotic second moment of the jump sigma_k^2")
# This is the expected squared change in position from state k.
# sigma_k^2 = (-2)^2*P_{k,k-2} + (-1)^2*P_{k,k-1} + (1)^2*P_{k,k+1} + (2)^2*P_{k,k+2}
# sigma_k^2 = 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k) + 4*(1/4)
# sigma_k^2 = 1 + 1/4 - c/k + 1/4 + c/k + 1 = 2.5
sigma_sq_B = 2.5
print(f"For large k, the second moment sigma_k^2 is constant, B = {sigma_sq_B}.")
print("-" * 40)

print("Step 3: Set up the transience condition")
# The criterion for transience is 2*A / B > 1.
# (2 * (2*c)) / 2.5 > 1
# 4*c / 2.5 > 1
# 8*c / 5 > 1
p_num = 8
p_den = 5
print(f"The condition for transience is (2*A)/B > 1, which becomes:")
print(f"({p_num}/{p_den}) * c > 1")
print("-" * 40)

print("Step 4: Solve for c and find the infimum")
# The inequality is 8/5 * c > 1  ==>  8*c > 5  ==> c > 5/8
print(f"Solving the inequality gives c > {p_den}/{p_num}.")
infimum = p_den / p_num
print(f"The set of c for which the chain is transient is ({infimum}, infinity).")
print(f"The infimum (greatest lower bound) of this set is {infimum}.")
print("-" * 40)

print("The final boundary equation is 8 * c = 5.")
print("As requested, the numbers in this final equation are:")
print(p_num)
print(p_den)

final_answer = p_den / p_num
print(f"\nThe final answer is {final_answer}.")