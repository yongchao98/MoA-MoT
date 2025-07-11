import sys

def solve():
    """
    Determines the maximal entropy H(x,y,z,s1,s2) subject to the given constraints.
    The method involves simplifying the entropy expression based on the constraints
    and then finding a configuration of random variables that achieves the maximum possible value.
    """

    print("Step 1: Simplify the joint entropy expression using the constraints.")
    print("=" * 70)
    print("The constraints of the form H(A | B, C) = 0 imply that A is a deterministic function of B and C.")
    print("From the given constraints, we can infer the following functional relationships:")
    print("1. H(z | s1, s2) = 0  =>  z is a function of (s1, s2).")
    print("2. H(x | s1, y) = 0  =>  x is a function of (s1, y).")
    print("3. H(y | x, s2) = 0  =>  y is a function of (x, s2).")
    print("\nBy substituting (3) into (2), we get that x depends on (s1, x, s2), which implies x is a function of (s1, s2).")
    print("Similarly, by substituting x into (3), we find that y is also a function of (s1, s2).")
    print("Therefore, all three variables x, y, and z are determined by s1 and s2.")
    print("\nBecause (x, y, z) are functionally dependent on (s1, s2), the entropy of the entire system of five variables")
    print("is equal to the entropy of the two base variables, s1 and s2:")
    print("H(x, y, z, s1, s2) = H(s1, s2)\n")

    print("Step 2: Establish an upper bound for H(s1, s2).")
    print("=" * 70)
    print("The chain rule for entropy states that H(s1, s2) = H(s1) + H(s2 | s1).")
    print("Since conditioning cannot increase entropy, H(s2 | s1) <= H(s2).")
    print("Therefore, H(s1, s2) <= H(s1) + H(s2).")
    print("The problem provides the constraints H(s1) <= 1 and H(s2) <= 1.")
    print("This gives us an upper bound for the entropy we want to maximize:")
    print("H(s1, s2) <= 1 + 1 = 2.\n")

    print("Step 3: Show that this upper bound of 2 is achievable.")
    print("=" * 70)
    print("To show the maximum is 2, we need to find a specific assignment for the random variables")
    print("that satisfies all constraints and achieves H(s1, s2) = 2.")
    print("H(s1, s2) = 2 is achieved if s1 and s2 are independent and H(s1) = 1, H(s2) = 1.")
    print("Let s1 and s2 be independent fair coin flips (i.e., i.i.d. Bernoulli(0.5) variables).")
    print("\nLet's propose a construction:")
    print("  s1, s2 ~ i.i.d. Bernoulli(0.5)")
    print("  x = s1")
    print("  y = s2")
    print("  z = s1 XOR s2\n")

    print("Checking if this construction satisfies all constraints:")
    print("  - Individual entropies: H(x)=H(s1)=1, H(y)=H(s2)=1, H(z)=1 (XOR of two iid bits is uniform). All <= 1. (OK)")
    print("  - H(s1|z,x) = H(s1|s1^s2, s1) = 0 since we are conditioning on s1. (OK)")
    print("  - H(s2|y,z) = H(s2|s2, s1^s2) = 0 since we are conditioning on s2. (OK)")
    print("  - H(x|s1,y) = H(s1|s1,s2)   = 0 since we are conditioning on s1. (OK)")
    print("  - H(y|x,s2) = H(s2|s1,s2)   = 0 since we are conditioning on s2. (OK)")
    print("  - H(z|s2,s1) = H(s1^s2|s1,s2) = 0 since z is a function of s1 and s2. (OK)")
    print("All constraints are satisfied.\n")

    print("Step 4: Final calculation.")
    print("=" * 70)
    print("With the proposed construction, the joint entropy is:")
    print("H(x, y, z, s1, s2) = H(s1, s2)")
    print("Since s1 and s2 are independent in this construction:")
    print("H(s1, s2) = H(s1) + H(s2)")
    h_s1 = 1
    h_s2 = 1
    max_entropy = h_s1 + h_s2
    print(f"The maximal entropy = {h_s1} + {h_s2} = {max_entropy}")

solve()
<<<2>>>