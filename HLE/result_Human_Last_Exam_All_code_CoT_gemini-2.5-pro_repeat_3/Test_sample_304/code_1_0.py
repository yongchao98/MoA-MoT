import math

def H_binary(p):
    """Calculates the entropy of a binary random variable."""
    if p == 0 or p == 1:
        return 0
    return -p * math.log2(p) - (1-p) * math.log2(1-p)

def solve_entropy_maximization():
    """
    Solves the entropy maximization problem by symbolic deduction.
    The steps of the deduction are printed to explain the reasoning.
    """
    
    print("### Step 1: Analyze the Constraints as Functional Dependencies ###")
    print("A constraint of the form H(A|B) = 0 means that A is a function of B.")
    print("From the problem statement, we have the following functional dependencies:")
    print("  (6) s1 is a function of (z, x)")
    print("  (7) s2 is a function of (y, z)")
    print("  (8) x is a function of (s1, y)")
    print("  (9) y is a function of (x, s2)")
    print("  (10) z is a function of (s1, s2)")
    print("-" * 50)

    print("### Step 2: Simplify the Objective Function ###")
    print("The goal is to maximize H(x, y, z, s1, s2).")
    print("From the dependencies, we can show that (x, y, z) are functions of (s1, s2):")
    print("  - From (8) and (9), x = f(s1, y) and y = g(x, s2).")
    print("  - Substituting y, we get x = f(s1, g(x, s2)), which means x is determined by (s1, s2). So, H(x|s1, s2) = 0.")
    print("  - Similarly, y is determined by (s1, s2). So, H(y|s1, s2) = 0.")
    print("  - From (10), z is determined by (s1, s2). So, H(z|s1, s2) = 0.")
    print("This means the entire vector (x, y, z) is a function of (s1, s2), so H(x, y, z | s1, s2) = 0.")
    print("\nUsing the chain rule of entropy:")
    print("  H(x, y, z, s1, s2) = H(s1, s2) + H(x, y, z | s1, s2)")
    print("  H(x, y, z, s1, s2) = H(s1, s2) + 0")
    print("Therefore, maximizing H(x, y, z, s1, s2) is equivalent to maximizing H(s1, s2).")
    print("-" * 50)
    
    print("### Step 3: Establish an Upper Bound ###")
    print("We know that entropy is sub-additive: H(s1, s2) <= H(s1) + H(s2).")
    print("The problem states that H(s1) <= 1 and H(s2) <= 1.")
    h_s1_max = 1
    h_s2_max = 1
    upper_bound = h_s1_max + h_s2_max
    print(f"So, H(s1, s2) <= {h_s1_max} + {h_s2_max} = {upper_bound}.")
    print(f"The maximal entropy is at most {upper_bound}.")
    print("-" * 50)
    
    print("### Step 4: Construct an Example to Achieve the Bound ###")
    print("To show that 2 is achievable, we need to find a set of random variables that satisfy all constraints.")
    print("Let s1 and s2 be independent, fair binary random variables (e.g., coin flips).")
    print(f"  - P(s1=0) = P(s1=1) = 0.5 => H(s1) = {H_binary(0.5)} = 1.")
    print(f"  - P(s2=0) = P(s2=1) = 0.5 => H(s2) = {H_binary(0.5)} = 1.")
    print("  - Since they are independent, H(s1, s2) = H(s1) + H(s2) = 1 + 1 = 2.")
    print("\nNow, define x, y, and z as functions of s1 and s2:")
    print("  Let x = s1")
    print("  Let y = s2")
    print("  Let z = s1 XOR s2")
    print("\nLet's verify all constraints with this construction:")
    print(f"  - H(x) = H(s1) = 1 <= 1 (OK)")
    print(f"  - H(y) = H(s2) = 1 <= 1 (OK)")
    print(f"  - z is also a fair coin flip, so H(z) = 1 <= 1 (OK)")
    print(f"  - H(s1) = 1 <= 1 and H(s2) = 1 <= 1 (OK)")
    print("  - H(s1|z,x) = H(s1|s1 XOR s2, s1) = 0, since x=s1 is given (OK)")
    print("  - H(s2|y,z) = H(s2|s2, s1 XOR s2) = 0, since y=s2 is given (OK)")
    print("  - H(x|s1,y) = H(s1|s1,s2) = 0, since s1 is given (OK)")
    print("  - H(y|x,s2) = H(s2|s1,s2) = 0, since s2 is given (OK)")
    print("  - H(z|s1,s2) = H(s1 XOR s2|s1,s2) = 0, since z is a function of (s1,s2) (OK)")
    print("\nAll constraints are satisfied.")
    print("-" * 50)

    print("### Step 5: Final Conclusion ###")
    print("The constructed distribution satisfies all constraints and achieves the upper bound.")
    print("The final calculation is:")
    h_s1 = 1
    h_s2 = 1
    max_H = h_s1 + h_s2
    print(f"max H(x,y,z,s1,s2) = H(s1) + H(s2) = {h_s1} + {h_s2} = {max_H}")
    return max_H

if __name__ == '__main__':
    final_answer = solve_entropy_maximization()
