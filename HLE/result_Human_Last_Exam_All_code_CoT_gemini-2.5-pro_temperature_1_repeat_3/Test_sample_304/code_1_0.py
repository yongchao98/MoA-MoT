def solve_entropy_maximization():
    """
    Solves the entropy maximization problem by simplifying the expression
    for the joint entropy and then finding its maximum value under the
    given constraints.
    """
    print("### Step-by-step Derivation ###\n")

    print("Step 1: Simplify the objective function H(x,y,z,s1,s2)")
    print("We start with the chain rule for entropy:")
    print("H(x,y,z,s1,s2) = H(x,y,z) + H(s1|x,y,z) + H(s2|x,y,z,s1)")
    print("\nThe constraint H(A|B) = 0 means A is fully determined by B.")
    print("From the constraint H(s1 | z,x) = 0, we know s1 is a function of (z,x).")
    print("Therefore, H(s1|x,y,z) <= H(s1|x,z) = 0, which implies H(s1|x,y,z) = 0.")
    print("\nSimilarly, from the constraint H(s2 | y,z) = 0, s2 is a function of (y,z).")
    print("Therefore, H(s2|x,y,z,s1) <= H(s2|y,z) = 0, which implies H(s2|x,y,z,s1) = 0.")
    print("\nSubstituting these results back, the joint entropy simplifies to:")
    print("H(x,y,z,s1,s2) = H(x,y,z) + 0 + 0 = H(x,y,z)\n")

    print("Step 2: Simplify H(x,y,z) using the dependency cycle")
    print("We can expand H(x,y,z) as: H(x,y,z) = H(x,y) + H(z|x,y)")
    print("Let's analyze the constraints involving z:")
    print("1. H(z | s2,s1) = 0  => z is a function of (s1, s2).")
    print("2. H(s1 | z,x) = 0  => s1 is a function of (z, x).")
    print("3. H(s2 | y,z) = 0  => s2 is a function of (y, z).")
    print("\nThese constraints imply z = f(s1, s2) = f(g(x,z), h(y,z)).")
    print("This means for any given pair of (x,y), z is determined by an equation involving only itself.")
    print("This implies z must be a deterministic function of (x,y), so H(z|x,y) = 0.")
    print("\nTherefore, the joint entropy further simplifies to:")
    print("H(x,y,z,s1,s2) = H(x,y,z) = H(x,y) + H(z|x,y) = H(x,y)\n")

    print("Step 3: Maximize H(x,y) to find the upper bound")
    print("The task is now to maximize H(x,y).")
    print("Using the chain rule and the property that conditioning does not increase entropy:")
    print("H(x,y) = H(x) + H(y|x) <= H(x) + H(y)")
    print("\nFrom the problem constraints, we have H(x) <= 1 and H(y) <= 1.")
    print("This gives us the final upper bound for the joint entropy.\n")

    print("### Final Calculation ###\n")
    # Define the variables for the final equation as per the constraints
    H_x_max = 1
    H_y_max = 1
    max_entropy = H_x_max + H_y_max
    
    print("The maximal entropy is bounded as follows:")
    print(f"H(x,y,z,s_1,s_2) = H(x,y)")
    print(f"               <= H(x) + H(y)")
    print(f"               <= {H_x_max} + {H_y_max}")
    print(f"               = {max_entropy}")

    print("\nThis maximum is achievable. For instance, let x and y be independent coin flips,")
    print("z = x XOR y, s1 = x, and s2 = y. This configuration satisfies all constraints")
    print(f"and results in a total entropy of {max_entropy}.")

if __name__ == "__main__":
    solve_entropy_maximization()