import math

# The problem is to find K, where the slice rank of a tensor is given by the formula (3/2^K)^n * e^(o(n)).
# Our analysis shows that the slice rank of the given tensor is exactly 2^n.
# By equating the asymptotic formula with our derived slice rank, we find the value for K.
# 2^n = (3 / (2^K))^n * e^(o(n))
# Taking the n-th root of both sides gives:
# 2 = (3 / (2^K)) * (e^(o(n)))^(1/n)
# As n -> infinity, the term (e^(o(n)))^(1/n) = e^(o(n)/n) -> e^0 = 1.
# So we have the equation:
# 2 = 3 / (2^K)
# 2 * 2^K = 3
# 2^(K+1) = 3
# K + 1 = log2(3)
# K = log2(3) - 1
# This is equivalent to K = log2(3/2).

K = math.log2(3.0) - 1.0
print(K)