import itertools

def compose(op1, op2):
    """
    Composes two affine transformations, where op1 is applied first.
    op1: (a1, b1) representing f1(x) = a1*x + b1
    op2: (a2, b2) representing f2(x) = a2*x + b2
    Returns a new transformation op_new = op2 o op1, where
    f_new(x) = f2(f1(x)) = a2*(a1*x + b1) + b2 = (a1*a2)*x + (a2*b1 + b2)
    """
    a1, b1 = op1
    a2, b2 = op2
    return (a1 * a2, a2 * b1 + b2)

# --- Horizontal Transformations ---
# We represent the transformations symbolically.
# R1(x) = -x is represented by (-1, 0)
# R3(x) = 2W - x is represented by (-1, 1) where 1 stands for 2W
r1 = (-1, 0)
r3 = (-1, 1)
h_ops_list = [r1, r1, r3]

# Get unique permutations of the horizontal reflection operators
h_perms = set(itertools.permutations(h_ops_list))

# Calculate the final transformation for each unique permutation
final_h_ops = set()
for perm in h_perms:
    # The sequence of reflections is perm[0], perm[1], ...
    # The transformation is T = T_n o ... o T_1
    final_op = perm[0]
    for i in range(1, len(perm)):
        final_op = compose(final_op, perm[i])
    final_h_ops.add(final_op)

num_h_transforms = len(final_h_ops)

# --- Vertical Transformations ---
# R2(y) = -y is represented by (-1, 0)
# R4(y) = 2H - y is represented by (-1, 2) where 2 stands for 2H
r2 = (-1, 0)
r4 = (-1, 2)
v_ops_list = [r2, r4]

# Get unique permutations of the vertical reflection operators
v_perms = set(itertools.permutations(v_ops_list))

# Calculate the final transformation for each unique permutation
final_v_ops = set()
for perm in v_perms:
    final_op = perm[0]
    for i in range(1, len(perm)):
        final_op = compose(final_op, perm[i])
    final_v_ops.add(final_op)

num_v_transforms = len(final_v_ops)

# The total number of ways is the product of the number of unique
# horizontal and vertical transformations.
total_ways = num_h_transforms * num_v_transforms

print("The number of unique paths is the product of the number of unique horizontal and vertical transformations.")
print(f"Number of unique horizontal transformations: {num_h_transforms}")
print(f"Number of unique vertical transformations: {num_v_transforms}")
print("Final calculation:")
print(f"{num_h_transforms} * {num_v_transforms} = {total_ways}")
