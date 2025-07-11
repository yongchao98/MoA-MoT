# Here is a code for Python 3.13.1.
m = True

# For set 'a':
# <A> is '^'. The initial list is [m ^ m, m], which evaluates to [1 ^ 1, 1] -> [0, 1].
# The loop f.append(f[-1] + f[-2]) generates a Fibonacci-like sequence.
# a becomes {0, 1, 2, 3, 5, 8, 13, 21, 34, 55}.
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([m ^ m, m]))

# For set 'b':
# <A> is '^', <C> is '*', <D> is '<<'.
# The initial list is [m ^ m] * (m << m) + [m] -> [0] * 2 + [1] -> [0, 0, 1].
# <B> is '|'. The index is ~(m | -~m) | m -> ~(1 | 2) | 1 -> ~3 | 1 -> -4 | 1 -> -3.
# The loop is f.append(f[-3] + f[-1] + f[-2]).
# b becomes {0, 1, 2, 4, 7, 13, 24, 44, 81, 149}.
b = set((lambda f: [f.append(f[~(m | -~m) | m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([m ^ m]*(m << m)+[m]))

# <F> is '&'. The intersection of a and b is {0, 1, 2, 13}.
final_set = b & a

# <E> is 'sum'. The sum of the intersection is 0 + 1 + 2 + 13 = 16.
result = sum(final_set)

# To meet the output requirement, we format the equation string.
sorted_list = sorted(list(final_set))
equation_str = ' + '.join(map(str, sorted_list))

print(f"{equation_str} = {result}")
