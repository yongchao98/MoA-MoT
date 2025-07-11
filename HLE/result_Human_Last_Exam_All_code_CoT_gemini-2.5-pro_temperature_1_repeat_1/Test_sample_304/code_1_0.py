import sympy

def solve_entropy_maximization():
    """
    Solves the entropy maximization problem by deriving an upper bound and
    providing a constructive example that achieves this bound.
    """

    # Represent the variables symbolically for clarity in the explanation
    H = sympy.Function('H')
    x, y, z, s1, s2 = sympy.symbols('x y z s_1 s_2')
    H_total = H(x, y, z, s1, s2)

    print("Problem: Determine the maximal entropy H(x,y,z,s1,s2) subject to the constraints:")
    print("H(x)<=1, H(y)<=1, H(z)<=1, H(s1)<=1, H(s2)<=1")
    print("H(s1|z,x)=0, H(s2|y,z)=0, H(x|s1,y)=0, H(y|x,s2)=0, H(z|s2,s1)=0")
    print("-" * 50)

    # Step 1: Derive the Upper Bound
    print("Step 1: Deriving an Upper Bound for H(x,y,z,s1,s2)\n")

    print("The condition H(A|B) = 0 means that A is a deterministic function of B.")
    print("From H(z|s2,s1)=0 and H(s1|z,x)=0, we can deduce that s1 is a function of (x,s2).")
    print("This implies that H(s1|x,s2) = 0.\n")

    print("We can expand the total entropy using the chain rule as follows:")
    print(f"{H_total} = H(x,s2) + H(y|x,s2) + H(s1|x,y,s2) + H(z|x,y,s1,s2)\n")

    print("Let's simplify this expression using the constraints:")
    print("1. H(y|x,s2) = 0, directly from the given constraints.")
    print("2. H(s1|x,y,s2) <= H(s1|x,s2). Since we derived H(s1|x,s2) = 0, this term is also 0.")
    print("3. H(z|x,y,s1,s2) <= H(z|s1,s2). Since H(z|s1,s2) = 0 is a given constraint, this term is also 0.\n")

    print("Therefore, the total entropy simplifies to:")
    print(f"{H_total} = H(x,s2)\n")

    print("Using the subadditivity property of entropy, H(A,B) <= H(A) + H(B):")
    print(f"H(x,s2) <= H(x) + H(s2)\n")

    print("From the constraints, we know H(x) <= 1 and H(s2) <= 1.")
    h_x_val = 1
    h_s2_val = 1
    upper_bound = h_x_val + h_s2_val
    print(f"So, {H_total} <= H(x) + H(s2) <= {h_x_val} + {h_s2_val} = {upper_bound}\n")
    print("This establishes an upper bound of 2 for the total entropy.")
    print("-" * 50)

    # Step 2: Construct an Example Achieving the Bound
    print("Step 2: Constructive Proof of Achievability\n")

    print("We need to show that the value 2 can be achieved. Consider the following construction:")
    print("Let s1 and x be independent random variables representing fair coin flips (Bernoulli(0.5)).")
    print("Then H(s1) = 1 and H(x) = 1.\n")

    print("Define the other variables as functions of s1 and x:")
    print("s2 = s1")
    print("z = s1")
    print("y = x XOR s1 (where XOR is the exclusive OR operation)\n")

    print("Let's verify this construction against the constraints:")
    print("- H(x)=1, H(s1)=1. Since s2=z=s1, H(s2)=H(z)=1. Since y=x XOR s1, H(y)=1. All <= 1 constraints are met.")
    print("- H(s1|z,x) = H(s1|s1,x) = 0. (OK)")
    print("- H(s2|y,z) = H(s1|y,s1) = 0. (OK)")
    print("- H(x|s1,y): Since y=x XOR s1, x can be determined by x=y XOR s1. So H(x|s1,y) = 0. (OK)")
    print("- H(y|x,s2) = H(y|x,s1): Since y=x XOR s1, y is determined by x and s1. So H(y|x,s1) = 0. (OK)")
    print("- H(z|s2,s1) = H(s1|s1,s1) = 0. (OK)\n")
    
    print("All constraints are satisfied.\n")
    
    print("Now, let's calculate the total entropy for this construction.")
    print("All variables are functions of the independent variables x and s1.")
    print(f"{H_total} = H(x,s1)")
    print("Since x and s1 are independent:")
    h_x_val_c = 1
    h_s1_val_c = 1
    final_entropy = h_x_val_c + h_s1_val_c
    print(f"{H_total} = H(x) + H(s1) = {h_x_val_c} + {h_s1_val_c} = {final_entropy}\n")

    # Step 3: Final Conclusion
    print("-" * 50)
    print("Conclusion:\n")
    print("We have shown that the maximal entropy is at most 2, and we have constructed a valid")
    print("case where the entropy is exactly 2.")
    print(f"Therefore, the maximal entropy is {final_entropy}.")

solve_entropy_maximization()
<<<2>>>