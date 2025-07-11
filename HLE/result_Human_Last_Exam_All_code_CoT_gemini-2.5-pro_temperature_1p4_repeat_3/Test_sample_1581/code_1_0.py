# The problem is to find the number of distinct homeomorphism classes
# for a compact connected metric space X, for which the n-th
# configuration space is disconnected for some n >= 2.

# Step 1: A compact connected metric space is called a continuum.
# Step 2: The n-th configuration space is F_n(X) = {(x_1, ..., x_n) in X^n | x_i != x_j for i != j}.
# Step 3: A topological result states that if F_n(X) is disconnected for some n >= 2,
# then F_2(X) is disconnected. So we only need to consider the case n=2.
# Step 4: A celebrated theorem by Karol Borsuk (1930) states that for any continuum X,
# the space F_2(X) is disconnected if and only if X is an arc.
# Step 5: An arc is a space homeomorphic to the closed interval [0, 1].
# Step 6: By definition, all arcs belong to the same homeomorphism class.
# Step 7: Therefore, any space X satisfying the condition must be homeomorphic to [0, 1].

# The number of distinct homeomorphism classes is the number of classes to which [0,1] belongs.
# That is just one class.
num_classes = 1

# Let's print the equation that represents the final answer.
# Let N be the number of homeomorphism classes.
print(f"N = {num_classes}")