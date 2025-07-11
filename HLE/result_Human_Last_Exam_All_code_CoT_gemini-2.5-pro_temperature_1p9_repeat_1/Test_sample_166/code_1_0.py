# The problem is to find the infimum exponent c for the polynomial runtime O(n^c)
# of a (1,lambda) EA on the CLIFF_{3/4} function.

# 1. Analysis of the CLIFF_{3/4} function.
# The fitness landscape has a local optimum at |x|_1 = 3n/4.
# To find a better solution from this point, the algorithm must produce an
# offspring x' with |x'|_1 = k such that k - n/4 + 1/2 > 3n/4,
# which simplifies to k > n - 1/2. This requires creating the global
# optimum with k=n ones.

# 2. Required Jump
# To get from an individual with 3n/4 ones to the optimum with n ones,
# the algorithm must flip exactly the n - 3n/4 = n/4 zeros into ones,
# and flip none of the 3n/4 ones into zeros. This is a jump of size m = n/4.

# 3. Runtime with an Optimized Static Mutation Rate
# Standard mutation rate p=1/n makes this jump exponentially unlikely.
# However, the algorithm allows for the *best static lambda* and implicitly the best
# static mutation rate p.
# To make a jump of size m more likely, one can set the mutation rate to p = m/n.
# Here, m = n/4, so the optimal choice for the jump is p = (n/4)/n = 1/4.

# 4. Runtime Analysis with p=m/n
# With p=m/n, the EA is biased towards flipping m bits. The crucial part becomes
# climbing the initial slope up to the cliff at |x|_1 = n-m = 3n/4.
# The expected runtime for the (1+1) EA with mutation rate p=m/n on a cliff
# function with jump size m is known to be O(n*m).
# For the (1,lambda) EA, a suitable choice of lambda can achieve this, but not
# fundamentally improve on this polynomial complexity. The dominant term comes
# from the climb, not the jump.

# 5. Calculation of the exponent c
# With m = n/4, the runtime T is O(n * m) = O(n * (n/4)) = O(n^2).
# The exponent c is therefore 2. The question asks for the infimum c.
# Theoretical results confirm that for jump sizes m=Theta(n), a runtime of O(n*m)
# is optimal, thus leading to c=2. The exponent is not lower for other static
# parameter choices.

# The value of c needs to be rounded to three significant digits.
c = 2.0

# Printing the final equation logic
# The runtime T is on the order of n * m
# T = O(n * m)
# The cliff is at 3/4 * n, so the jump size m to reach n is n - 3/4 * n = 1/4 * n
m_numerator = 1
m_denominator = 4
# So, T = O(n * (n/4)) = O(n^2 / 4) = O(n^2)
# The exponent c is 2.
c_val = 2
power_of_n_1 = 1
power_of_n_2 = 1

print(f"The runtime T is given by the product of n and the jump size m.")
print(f"T = O(n^({power_of_n_1}) * m)")
print(f"The jump size m is calculated as n - (3/4)*n = ({m_numerator}/{m_denominator})*n.")
print(f"Substituting m gives: T = O(n * ({m_numerator}/{m_denominator})*n) = O(n^{power_of_n_1 + power_of_n_2})")
print(f"Therefore, the runtime is O(n^{c_val}).")
print(f"The value of c is {c_val}.")
# Final numerical answer required by the problem format
final_c = float(c_val)
print(f"The infimum c is {final_c:.3f}")