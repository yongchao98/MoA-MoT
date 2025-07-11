import numpy as np

def check_mixed_dimension_learnability(op_name, op_func):
    """
    For a mixed-dimension operator op(h1[i], h2[j]), its learnability requires it to be
    additively separable, meaning op(x, y) = f(x) + g(y) for some functions f and g.
    This function checks if such f and g can exist by solving a system of linear equations.

    The system of equations is derived from the four possible boolean inputs for (x, y):
    1. f(0) + g(0) = op(0, 0)
    2. f(1) + g(0) = op(1, 0)
    3. f(0) + g(1) = op(0, 1)
    4. f(1) + g(1) = op(1, 1)

    This can be written in matrix form A*v = b, where v = [f(0), f(1), g(0), g(1)].
    """
    
    # The coefficient matrix A is the same for all operators.
    A = np.array([
        # f(0), f(1), g(0), g(1)
        [1,    0,    1,    0],   # Equation 1: f(0) + g(0)
        [0,    1,    1,    0],   # Equation 2: f(1) + g(0)
        [1,    0,    0,    1],   # Equation 3: f(0) + g(1)
        [0,    1,    0,    1]    # Equation 4: f(1) + g(1)
    ])

    # The vector b contains the results of the operator for the four input pairs.
    b = np.array([
        op_func(0, 0),  # op(0, 0)
        op_func(1, 0),  # op(1, 0)
        op_func(0, 1),  # op(0, 1)
        op_func(1, 1)   # op(1, 1)
    ])

    print(f"--- Analyzing Operator: {op_name} ---")
    print(f"The required outputs for op(0,0), op(1,0), op(0,1), op(1,1) form the equation's right side: {b.tolist()}")

    # A system Ax=b has a solution if and only if rank(A) == rank([A|b]).
    # We can use numpy to check the ranks.
    rank_A = np.linalg.matrix_rank(A)
    # Construct the augmented matrix [A|b] to find its rank.
    augmented_A = np.c_[A, b]
    rank_augmented_A = np.linalg.matrix_rank(augmented_A)

    if rank_A == rank_augmented_A:
        is_learnable = True
        print("Result: System has a solution. The operator IS additively separable and thus theoretically learnable.")
    else:
        is_learnable = False
        print("Result: System is inconsistent (no solution). The operator is NOT additively separable and thus CANNOT be learned.")
    
    print("-" * 35 + "\n")
    return is_learnable

if __name__ == "__main__":
    # Define the mixed-dimension operators
    mixed_operators = {
        "X' (Mixed XOR)": lambda x, y: x ^ y,
        "C' (Mixed Conjunction)": lambda x, y: x & y,
        "D' (Mixed Disjunction)": lambda x, y: x | y,
        "E' (Mixed Equivalence)": lambda x, y: 1 - (x ^ y),
        "I' (Mixed Implication)": lambda x, y: (1 - x) | y,
    }

    not_learnable_list = []
    
    print("Element-wise operators (X, C, D, E, I) are all learnable because the heuristic")
    print("provides features like |h1-h2| (XOR) and h1*h2 (AND) directly.\n")
    
    print("Checking which mixed-dimension operators are NOT learnable:")
    print("=" * 60)

    for name, func in mixed_operators.items():
        op_code = name.split(" ")[0]
        if not check_mixed_dimension_learnability(name, func):
            not_learnable_list.append(op_code)

    print("=" * 60)
    print("Summary:")
    print("The operators that CANNOT be learned are those that are not additively separable.")
    print("Our analysis shows this applies to all mixed-dimension operators.")
    print(f"\nFinal list of non-learnable operators: {''.join(sorted(not_learnable_list))}")
    
    print("\n<<<I>>>")
