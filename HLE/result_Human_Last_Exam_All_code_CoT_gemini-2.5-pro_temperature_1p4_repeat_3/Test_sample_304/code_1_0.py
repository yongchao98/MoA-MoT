import math

def H(prob):
    """Calculates Shannon entropy for a binary variable with probability prob."""
    if prob == 0 or prob == 1:
        return 0
    return -prob * math.log2(prob) - (1 - prob) * math.log2(1 - prob)

def main():
    """
    Solves the entropy maximization problem.
    """
    print("Step 1: Simplify the objective function H(x,y,z,s1,s2).")
    print("We use the chain rule for entropy: H(A,B,C) = H(A) + H(B|A) + H(C|A,B).")
    print("The total entropy is H(x,y,z,s1,s2) = H(x,s2) + H(y,z,s1 | x,s2).")
    print("Let's expand the conditional part: H(y,z,s1 | x,s2) = H(y|x,s2) + H(z|x,s2,y) + H(s1|x,s2,y,z).")
    print("\nAnalyzing each term using the given constraints:")
    print("1. From the constraint H(y | x,s2) = 0, the first term is 0.")
    print("2. The property H(A|B,C) <= H(A|B) and the constraint H(z | s2,s1) = 0 imply H(z|x,s2,y) <= H(z|s2,s1) = 0. So the second term is 0.")
    print("3. Similarly, H(s1|x,s2,y,z) <= H(s1|z,x) = 0. So the third term is 0.")
    print("Therefore, H(y,z,s1 | x,s2) = 0 + 0 + 0 = 0.")
    print("This simplifies the total entropy to H(x,y,z,s1,s2) = H(x,s2).")

    print("\nStep 2: Find an upper bound for the simplified entropy H(x,s2).")
    print("Using the subadditivity property of entropy: H(x,s2) <= H(x) + H(s2).")
    print("From the given constraints, we have H(x) <= 1 and H(s2) <= 1.")
    H_x_val = 1
    H_s2_val = 1
    upper_bound = H_x_val + H_s2_val
    print(f"So, H(x,y,z,s1,s2) <= H(x) + H(s2) <= {H_x_val} + {H_s2_val} = {upper_bound}.")
    print(f"The maximal entropy is at most {upper_bound}.")

    print("\nStep 3: Construct an example to show this upper bound is achievable.")
    print("Let x and y be independent random variables from a fair coin flip (Bernoulli(0.5)).")
    print("Then H(x) = 1 and H(y) = 1.")
    print("Define the other variables as follows:")
    print("s1 = x")
    print("s2 = y")
    print("z = x XOR y")
    print("\nVerifying the constraints for this construction:")
    print(f"H(x) = {H(0.5)} <= 1. (Ok)")
    print(f"H(y) = {H(0.5)} <= 1. (Ok)")
    print(f"H(s1) = H(x) = {H(0.5)} <= 1. (Ok)")
    print(f"H(s2) = H(y) = {H(0.5)} <= 1. (Ok)")
    print(f"H(z) = H(x XOR y) = {H(0.5)} <= 1 (since x,y are iid Bernoulli(0.5), their XOR is also Bernoulli(0.5)). (Ok)")
    print("H(s1 | z,x) = H(x | x XOR y, x) = 0. (Ok, since x is given).")
    print("H(s2 | y,z) = H(y | y, x XOR y) = 0. (Ok, since y is given).")
    print("H(x | s1,y) = H(x | x,y) = 0. (Ok, since x is given).")
    print("H(y | x,s2) = H(y | x,y) = 0. (Ok, since y is given).")
    print("H(z | s2,s1) = H(x XOR y | y,x) = 0. (Ok, since z is a function of x and y).")
    print("All constraints are satisfied.")
    
    print("\nStep 4: Calculate the entropy for the constructed example.")
    print("H(x,y,z,s1,s2) = H(x,s2) = H(x,y).")
    print("Since x and y are independent in our construction:")
    final_h_x = H(0.5)
    final_h_y = H(0.5)
    final_entropy = final_h_x + final_h_y
    print(f"H(x,y) = H(x) + H(y) = {final_h_x} + {final_h_y} = {final_entropy}.")

    print("\nConclusion:")
    print("We have shown that the entropy is at most 2, and we have found a valid configuration where the entropy is exactly 2.")
    print(f"Therefore, the maximal entropy H(x,y,z,s1,s2) is {final_entropy}.")

if __name__ == '__main__':
    main()