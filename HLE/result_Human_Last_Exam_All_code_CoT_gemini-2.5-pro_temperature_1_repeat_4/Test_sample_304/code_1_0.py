import math

def solve():
    """
    Calculates the maximal entropy based on the derivation.

    The problem is to find the maximum of H(x, y, z, s1, s2) subject to a set of constraints.
    The constraints are:
    1. H(v) <= 1 for v in {x, y, z, s1, s2}
    2. H(s1|z,x)=0, H(s2|y,z)=0, H(x|s1,y)=0, H(y|x,s2)=0, H(z|s2,s1)=0

    Step-by-step derivation:
    1. Let H_total = H(x, y, z, s1, s2).
    2. Using the chain rule: H_total = H(s1, s2) + H(x, y, z | s1, s2).
    3. The constraint H(z|s2,s1)=0 means z is a function of (s1, s2).
       So, H(z|s1,s2)=0. The conditional entropy term becomes:
       H(x, y, z | s1, s2) = H(z|s1,s2) + H(x, y | s1, s2, z) = H(x, y | s1, s2, z).
    4. The remaining constraints imply that given (s1, s2, z), both x and y are uniquely determined.
       For instance, H(s1|z,x)=0 implies x is a function of (s1,z). Since s1 and z are given in the conditioning, x is determined.
       Similarly, H(s2|y,z)=0 implies y is a function of (s2,z), so y is also determined.
       Thus, the conditional entropy H(x, y | s1, s2, z) = 0.
    5. This simplifies the total entropy to H_total = H(s1, s2).
    6. To maximize H(s1, s2), we use the property H(s1, s2) <= H(s1) + H(s2).
    7. From the constraints H(s1) <= 1 and H(s2) <= 1, the maximum value is bounded by 1 + 1 = 2.
    8. This maximum is achievable when s1 and s2 are independent and have maximum entropy.
       For example, let s1 and s2 be independent fair coin flips.
       H(s1) = 1 bit.
       H(s2) = 1 bit.
    9. The maximal entropy is therefore H(s1) + H(s2).
    """

    # Maximum entropy for s1
    H_s1_max = 1

    # Maximum entropy for s2
    H_s2_max = 1

    # The maximal joint entropy H(x,y,z,s1,s2) is equal to H(s1,s2)
    # which is maximized when s1 and s2 are independent.
    # H_max = max(H(s1)) + max(H(s2))
    max_entropy = H_s1_max + H_s2_max

    print(f"The derivation shows the maximal entropy is H(s1) + H(s2).")
    print(f"Given H(s1) <= {H_s1_max} and H(s2) <= {H_s2_max}, the maximal value is achieved when H(s1) and H(s2) are at their maximums and s1 and s2 are independent.")
    print(f"Maximal Entropy = H(s1) + H(s2) = {H_s1_max} + {H_s2_max} = {max_entropy}")

solve()