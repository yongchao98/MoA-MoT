import math

# The problem is to find the Lebesgue measure of the set S of initial points x_0
# for which the sequence x_{n+1} = f(x_n) has exactly 7 distinct values.
# The function is f(x) = (2*x + sin(2*pi*x))/3.

# 1. A sequence with a finite number of distinct values must be eventually periodic.
#    Let the number of distinct values be N. This means the sequence has a tail of
#    some length k >= 0 and then enters a cycle of some length p >= 1.
#    The total number of distinct values is N = k + p.
#    In our case, N=7.

# 2. The set S is therefore the set of all points x_0 that are eventually periodic,
#    with k+p=7. This includes purely periodic points (k=0, p=7) and
#    points that are strictly pre-periodic (k>0, p>=1).

# 3. For a C^1 function on a compact interval, like our f(x) on [0,1], the set
#    of periodic points is countable. This is because for any period p, the set
#    of points satisfying f^p(x) = x is the set of roots of a function g(x) = f^p(x) - x.
#    Unless g(x) is identically zero on an interval, the set of roots is a set of
#    isolated points, which must be finite on a compact interval. Our f(x) is real
#    analytic and not the identity function, so f^p(x) is not x on any interval.
#    Thus, the set of all periodic points is a countable union of finite sets,
#    which is a countable set.

# 4. The set S consists of periodic points and their pre-images.
#    The function f(x) has a finite number of pre-images for any point y (at most 3).
#    The set of pre-images of a countable set is itself countable.
#    Therefore, the set S is a countable set.

# 5. A fundamental property of the Lebesgue measure is that any countable set of
#    real numbers has a measure of zero.

# 6. Therefore, the Lebesgue measure of S is 0.
lebesgue_measure_of_S = 0

# The question asks for this measure multiplied by 10^6.
result = lebesgue_measure_of_S * 10**6

print(int(result))