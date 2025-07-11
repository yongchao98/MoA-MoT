import math

def entropy(p):
    """Calculates entropy for a binary variable with probability p."""
    if p == 0 or p == 1:
        return 0
    return -p * math.log2(p) - (1-p) * math.log2(1-p)

def main():
    """
    Solves for the maximal entropy based on the provided constraints.
    """

    # 1. Simplify the objective function H(x,y,z,s1,s2)
    # The joint entropy can be expanded using the chain rule:
    # H(x,y,z,s1,s2) = H(s1,s2) + H(z|s1,s2) + H(x,y|s1,s2,z)
    # The constraint H(z|s1,s2) = 0 simplifies this to:
    # H(x,y,z,s1,s2) = H(s1,s2) + H(x,y|s1,s2)
    # Expanding H(x,y|s1,s2) = H(x|s1,s2) + H(y|x,s1,s2)
    # The constraint H(y|x,s2) = 0 implies H(y|x,s1,s2) = 0.
    # So, the final simplified expression for the total entropy is:
    # H(x,y,z,s1,s2) = H(s1,s2) + H(x|s1,s2)

    # 2. Construct a set of variables to satisfy the constraints and maximize entropy.
    # To maximize H(s1,s2) + H(x|s1,s2), we should try to maximize each part.
    # H(s1,s2) <= H(s1) + H(s2) <= 1 + 1 = 2.
    # We choose s1 and s2 to be independent, fair binary random variables (coin flips).
    
    # H(s1) for a fair coin (p=0.5)
    H_s1 = entropy(0.5)
    
    # H(s2) for a fair coin (p=0.5)
    H_s2 = entropy(0.5)
    
    # Since s1 and s2 are independent, H(s2|s1) = H(s2)
    H_s2_g_s1 = H_s2
    
    # H(s1, s2) = H(s1) + H(s2|s1) = 1 + 1 = 2
    H_s1s2 = H_s1 + H_s2_g_s1
    
    # Now, we define the other variables to satisfy the constraints. Let:
    # x = s1
    # y = s2
    # z = s1 XOR s2
    # This construction satisfies all problem constraints as shown in the explanation.

    # 3. Calculate the remaining term H(x|s1,s2).
    # Since x is defined as s1, knowing s1 and s2 means x is fully determined.
    # Therefore, the conditional entropy H(x|s1,s2) is 0.
    H_x_g_s1s2 = 0
    
    # 4. Calculate the total entropy.
    # H_total = H(s1,s2) + H(x|s1,s2)
    result = H_s1s2 + H_x_g_s1s2

    # We can also use an alternative decomposition of the objective function:
    # H(x,y,z,s1,s2) = H(s1) + H(s2|s1) + H(x|s1,s2) + H(y|x,s1,s2) + H(z|x,y,s1,s2)
    # The constraints H(y|x,s2)=0 and H(z|s1,s2)=0 make the last two terms zero.
    # So, H_total = H(s1) + H(s2|s1) + H(x|s1,s2)
    
    print("To find the maximal entropy, we can decompose the joint entropy and use a construction that satisfies the constraints.")
    print("A valid construction is to let s1 and s2 be independent fair coin flips, and define x=s1, y=s2, and z=s1 XOR s2.")
    print("For this construction, all constraints are met.")
    print("\nThe joint entropy H(x,y,z,s1,s2) can be calculated as H(s1) + H(s2|s1) + H(x|s1,s2).")
    print(f"H(s1) = {H_s1}")
    print(f"H(s2|s1) = H(s2) = {H_s2_g_s1} (due to independence)")
    print(f"H(x|s1,s2) = H(s1|s1,s2) = {H_x_g_s1s2}")
    print("\nThe final equation is:")
    print(f"H(x,y,z,s1,s2) = H(s1) + H(s2|s1) + H(x|s1,s2) = {H_s1} + {H_s2_g_s1} + {H_x_g_s1s2} = {result}")

if __name__ == "__main__":
    main()
