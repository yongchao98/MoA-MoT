import sys

def solve():
    """
    Solves the theoretical group theory problem and prints the answer.

    The reasoning is as follows:
    (a) The group hat(G) is the p-localization of G. For the class of solvable groups to which G belongs,
        the p-localization exists and is unique up to isomorphism due to its universal property.
        Therefore, a minimal such group exists and is unique.

    (b) Let dl(X) be the derived length of group X.
        - The existence of a solvable series of length n for G implies dl(G) <= n.
        - Localization has the property dl(hat(G)) <= dl(G). Thus, dl(hat(G)) <= n.
        - To show n is the maximum, we need to show it's achievable. This happens if we can find a group G
          satisfying the conditions for which dl(G) = n and dl(hat(G)) = dl(G).
        - The equality dl(hat(G)) = dl(G) holds if G is p'-torsion-free.
        - G is p'-torsion-free if all its composition factors B_i are. The problem requires this for i < n.
          We are free to construct a G where B_n is also p'-torsion-free.
        - An iterated wreath product of a p'-torsion-free abelian group (like Z) gives a group G
          with dl(G)=n that satisfies all conditions.
        - For this group, dl(hat(G)) = n.
        - Thus, the maximum possible derived length is n.
    """
    
    # Part (a): Existence and uniqueness of the minimal group
    answer_a = "Yes"
    
    # Part (b): Maximum possible derived length of the group.
    # The final equation is dl_max = n
    answer_b_expression = "n"
    
    print(f"(a) [{answer_a}]; (b) [{answer_b_expression}]")

solve()
