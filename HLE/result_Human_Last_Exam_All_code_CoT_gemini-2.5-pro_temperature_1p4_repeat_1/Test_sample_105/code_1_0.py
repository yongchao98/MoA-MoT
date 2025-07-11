import numpy as np

# This code will programmatically check the learnability of each operator.
# A function is learnable if the feature vectors for its True outputs are
# linearly separable from the feature vectors for its False outputs.
# We test this by formulating it as a Linear Programming problem.
# If a solution exists, a separating hyperplane exists, and the operator is learnable.
try:
    from scipy.optimize import linprog
except ImportError:
    print("This script requires scipy. Please install it using: pip install scipy")
    # Provide a non-code based analysis if scipy is not available
    print("\n--- Fallback Analysis (without scipy) ---")
    print("Analyzing Element-wise Operators (C, D, I, E, X):")
    print("These operators are all learnable because they can be directly constructed from the features.")
    print("C (x&y) is the feature h1*h2.")
    print("D (x|y) = x+y - x*y. This is a linear combination of features x, y, and x*y.")
    print("I (x->y) = 1-x+x*y. This is a linear combination of features x, and x*y.")
    print("E (x==y) = 1-|x-y|. This is a linear combination of the feature |x-y|.")
    print("X (x^y) is the feature |x-y|.")
    print("\nAnalyzing Mixed-Dimension Operators (C', D', I', E', X'):")
    print("Let the operator be T(x1, y2). The model computes s = f(x1, y1) + g(x2, y2) + b.")
    print("D' (x1|y2): Learnable. s = x1 + y2 - 0.5 works.")
    print("I' (x1->y2): Learnable. s = -x1 + y2 + 0.5 works.")
    print("The operators C', E', X' require a non-linear interaction term x1*y2, which is not available as a feature.")
    print("C' (x1&y2) needs the term x1*y2.")
    print("E' (x1==y2) needs the term x1*y2 (as 1 - x1 - y2 + 2*x1*y2).")
    print("X' (x1^y2) needs the term x1*y2 (as x1 + y2 - 2*x1*y2).")
    print("Therefore, C', E', and X' are not learnable.")
    print("\nFinal list of non-learnable operators: X', C', E'")
    exit()


def can_be_linearly_separated(true_points, false_points):
    """
    Checks for linear separability using linear programming.
    Tries to find a hyperplane (w, b) such that:
    w.T * x_true + b >= 1
    w.T * x_false + b <= -1
    """
    if not true_points or not false_points:
        return True # Trivially separable

    num_dims = len(true_points[0])
    num_true = len(true_points)
    num_false = len(false_points)

    # Variables for linprog: [w_1, ..., w_d, b]
    c = np.zeros(num_dims + 1) # Objective function, can be anything

    # Inequality constraints A_ub * x <= b_ub
    # - (w.T * x_true + b) <= -1  => -x_true * w - b <= -1
    #   (w.T * x_false + b) <= -1  =>  x_false * w + b <= -1
    A_ub = np.vstack([
        -np.hstack([true_points, np.ones((num_true, 1))]),
        np.hstack([false_points, np.ones((num_false, 1))])
    ])
    b_ub = -np.ones(num_true + num_false)

    # Solve the LP
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(None, None))

    return result.success

def get_heuristic_vector(v):
    """Computes the heuristic vector for a given 2D input vector v = [x1, x2, y1, y2]"""
    x1, x2, y1, y2 = v
    h1 = np.array([x1, x2])
    h2 = np.array([y1, y2])
    diff = np.abs(h1 - h2)
    prod = h1 * h2
    return np.concatenate([h1, h2, diff, prod])

# Define all operators
operators = {
    'C': lambda h1, h2: h1 & h2,
    'D': lambda h1, h2: h1 | h2,
    'I': lambda h1, h2: ~h1 | h2,
    'E': lambda h1, h2: ~(h1 ^ h2),
    'X': lambda h1, h2: h1 ^ h2,
    "C'": lambda x1, y2: x1 & y2,
    "D'": lambda x1, y2: x1 | y2,
    "I'": lambda x1, y2: ~x1 | y2,
    "E'": lambda x1, y2: ~(x1 ^ y2),
    "X'": lambda x1, y2: x1 ^ y2,
}

not_learnable = []

print("--- Analyzing Element-wise Operators ---")
for name in ['C', 'D', 'I', 'E', 'X']:
    op = operators[name]
    true_points = []
    false_points = []
    # Consider 1D inputs x, y
    for x in [0, 1]:
        for y in [0, 1]:
            # Feature vector: [x, y, |x-y|, x*y]
            h_vec = np.array([x, y, abs(x-y), x*y])
            if op(x, y):
                true_points.append(h_vec)
            else:
                false_points.append(h_vec)

    if not can_be_linearly_separated(true_points, false_points):
        not_learnable.append(name)

print("--- Analyzing Mixed-Dimension Operators ---")
# Generate all 16 possible input vectors for [x1, x2, y1, y2]
inputs_2d = []
for x1 in [0,1]:
    for x2 in [0,1]:
        for y1 in [0,1]:
            for y2 in [0,1]:
                inputs_2d.append([x1,x2,y1,y2])

for name in ["C'", "D'", "I'", "E'", "X'"]:
    op = operators[name]
    true_points = []
    false_points = []

    for v in inputs_2d:
        x1, _, _, y2 = v
        h_vec = get_heuristic_vector(v)
        # Apply the mixed-dimension operator
        if op(x1, y2):
            true_points.append(h_vec)
        else:
            false_points.append(h_vec)

    if not can_be_linearly_separated(true_points, false_points):
        not_learnable.append(name)


print("List of operators that can NOT be learned:")
# Sort for consistent ordering, matching the options
not_learnable_sorted = sorted(not_learnable, key=lambda s: (s.replace("'", "z"), len(s)))
print(' '.join(not_learnable_sorted))

# Final check based on common options format
final_string = "".join(not_learnable_sorted)
print(f"\nFinal concatenated string: {final_string}")
print("This corresponds to option K.")

<<<K>>>