# The problem asks for the smallest possible number of isometric embeddings
# of a finite ultrametric space X into a Banach space B.

# Let's denote an isometric embedding as a function f: X -> B such that
# ||f(x) - f(y)|| = d(x, y) for all x, y in X.

# If f is an isometric embedding, let's consider the function g(x) = -f(x).
# The distance preservation for g is shown by:
# ||g(x) - g(y)|| = ||-f(x) - (-f(y))|| = ||-(f(x) - f(y))|| = ||f(x) - f(y)|| = d(x, y).
# So, g is also an isometric embedding.

# The embedding g is distinct from f, unless f(x) = -f(x) for all x, which
# would mean f(x) = 0. This "zero embedding" is only possible if all distances
# in X are 0, i.e., X is a trivial single-point space.

# For any non-trivial finite ultrametric space, an embedding f exists into some
# Banach space, and it won't be the zero map.
# Therefore, if one such embedding (f) exists, its distinct negative partner (-f) must also exist.
# This gives us a count of at least two embeddings.

# One embedding is the original one found.
first_embedding = 1

# The second is its negative, which is also a valid and distinct embedding.
second_embedding = 1

# The smallest possible number of embeddings is the sum of these two.
total = first_embedding + second_embedding

print(f"{first_embedding} + {second_embedding} = {total}")