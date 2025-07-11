import numpy as np

def h(p):
    """Calculates the entropy of a Bernoulli variable."""
    if p == 0 or p == 1:
        return 0
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

def main():
    """
    This function demonstrates the calculation for the maximal entropy based on the derived result.
    The derivation shows that the maximal entropy is H(x) + H(y) under a valid construction
    that sets H(z|x,y) = 0.
    """
    
    # We found a construction that satisfies all constraints:
    # x ~ Bernoulli(0.5), y ~ Bernoulli(0.5), x and y are independent.
    # s1 = x, s2 = y, z = x XOR y.
    
    # Entropies of the independent variables
    # H(x) for p=0.5 is 1
    # H(y) for p=0.5 is 1
    h_x = h(0.5)
    h_y = h(0.5)
    
    # From our derivation, the maximal entropy H(x,y,z,s1,s2) = H(x,y)
    # Since x and y are independent, H(x,y) = H(x) + H(y)
    max_entropy = h_x + h_y
    
    print(f"The derivation leads to the inequality H(x,y) + 0.5 * H(z|x,y) <= 2.")
    print(f"The constraints force H(z|x,y) = 0.")
    print(f"This makes the problem equivalent to maximizing H(x,y) subject to H(x)<=1 and H(y)<=1.")
    print(f"The maximum value of H(x,y) is H(x) + H(y).")
    print(f"With H(x) <= 1 and H(y) <= 1, the maximum is achieved when H(x)=1 and H(y)=1.")
    print(f"Final Equation: H_max = H(x) + H(y)")
    print(f"H_max = {h_x} + {h_y} = {max_entropy}")
    print(f"The maximal entropy is {max_entropy}")

if __name__ == "__main__":
    main()
