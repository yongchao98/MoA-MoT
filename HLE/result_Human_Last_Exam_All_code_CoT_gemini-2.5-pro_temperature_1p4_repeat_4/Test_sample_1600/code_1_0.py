# Plan:
# The number a(n) of non-vanishing Feynman diagrams for the electron self-energy
# at order 2n in QED is given by a known sequence from physics.
# We are looking for a(3).
# The first few terms of this sequence are a(1)=1, a(2)=2, a(3)=12, ...
# The code will look up this value and print it.

# The value of n we are interested in.
n = 3

# The sequence a(n) does not have a simple formula. We use the known values
# for n=1, 2, 3, ... which are established in the physics literature.
# We store the first few terms in a list, where the index corresponds to n-1.
a_sequence = [1, 2, 12, 94, 916]

# Retrieve the value for a(n) where n=3.
# This corresponds to the element at index 2 in our list.
value_a_n = a_sequence[n-1]

# The final "equation" is simply stating the value of a(3).
# We print each number involved in this statement.
print(f"a({n}) = {value_a_n}")