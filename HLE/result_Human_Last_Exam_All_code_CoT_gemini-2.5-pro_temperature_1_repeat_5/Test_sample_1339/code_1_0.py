import sys

def solve_group_theory_problem():
    """
    Solves the given group theory problem regarding derived length.

    The reasoning is as follows:

    Part (a): Does there exist a unique minimal group \hat{G}?
    The conditions on \hat{G} (that it contains G and every p-nonsingular system of equations over G has a solution in \hat{G})
    define the p-completion, or p-localization, of G. For the class of solvable groups described, the existence of such a
    completion is a standard result in group theory. This group \hat{G} is constructed to be minimal and is unique
    up to a unique isomorphism that is the identity on G. This uniqueness is guaranteed by the universal property that
    defines the construction. Therefore, the answer is "Yes".

    Part (b): What is the maximum possible derived length of \hat{G}?
    Let dl(H) denote the derived length of a group H.
    1.  The p-completion functor preserves the derived length for the class of solvable groups considered here.
        This means dl(\hat{G}) = dl(G). So, we need to find the maximum possible derived length of G.

    2.  The group G has a subnormal series G = G_1 \triangleleft G_2 \triangleleft \dots \triangleleft G_n \triangleleft G_{n+1} = {1}
        where each quotient G_i / G_{i+1} is abelian. This is a solvable series of length n.
        Let's relate this to the derived series G^(k).
        - G^(0) = G = G_1.
        - Since G_1/G_2 is abelian, the commutator subgroup [G_1, G_1] must be contained in G_2. So, G^(1) \subseteq G_2.
        - Inductively, if G^(k-1) \subseteq G_k, then G^(k) = [G^(k-1), G^(k-1)] \subseteq [G_k, G_k].
        - Since G_k/G_{k+1} is abelian, [G_k, G_k] \subseteq G_{k+1}.
        - Thus, G^(k) \subseteq G_{k+1}.
        - For k=n, this gives G^(n) \subseteq G_{n+1} = {1}. So, G^(n) = {1}.
        - This proves that the derived length of G is at most n, i.e., dl(G) <= n.

    3.  To find the maximum possible derived length, we must check if this upper bound `n` can be achieved.
        It is indeed possible to construct a group G with derived length exactly n that meets the criteria.
        For example, one can construct a solvable group G of derived length n whose derived factors G^(i-1)/G^(i)
        are p'-torsion-free for i < n (e.g., they can be made free abelian), and the last non-trivial term of the derived series G^(n-1) is a suitable abelian group.
        Such a group G has dl(G) = n, and its derived series G \triangleright G^(1) \triangleright \dots \triangleright G^(n-1) \triangleright {1}
        can be taken as the series G_1, G_2, ..., G_{n+1}, which satisfies all the conditions of the problem.

    4.  Since dl(G) <= n and dl(G) can be equal to n, the maximum possible derived length for G is n.
        Consequently, the maximum possible derived length for \hat{G} is also n.
    """

    # Answer for part (a)
    answer_a = "Yes"

    # Answer for part (b)
    # The derived length depends on the parameter 'n' from the problem description.
    answer_b = "n"

    print(f"(a) {answer_a}; (b) {answer_b}")

solve_group_theory_problem()