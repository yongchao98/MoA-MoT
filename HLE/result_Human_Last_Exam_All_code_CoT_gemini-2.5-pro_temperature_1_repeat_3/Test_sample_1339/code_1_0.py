import sys

def solve_group_theory_problem():
    """
    Solves the given group theory problem and prints the answer.

    The problem asks two questions about a group `hat(G)` derived from a group G.
    G has a subnormal series G = G_1 > G_2 > ... > G_n > G_{n+1} = {1}
    with abelian factors B_i = G_i / G_{i+1}. The factors B_i for i < n are p'-torsion-free.
    `hat(G)` is a group containing G where any p-nonsingular system over G has a solution.

    The variable 'n' is symbolic and comes from the problem description.
    """

    # Part (a): Does there exist a unique minimal group `hat(G}`?

    # The group `hat(G)` described is the p-localization of G. For any group G,
    # a p-localization exists and is unique up to isomorphism.
    # The group G given in the problem is poly-abelian, which means it is solvable.
    # For a solvable group G, the localization map l: G -> hat(G) is injective.
    # This means we can identify G with its image in `hat(G)`, so `hat(G)` contains G.
    # The universal property of p-localization ensures that `hat(G)` is the minimal
    # such group. Any other group with the required property must contain an
    # isomorphic copy of `hat(G)`.
    # Therefore, the answer is Yes.
    answer_a = "Yes"

    # Part (b): What is the maximum possible derived length of `hat(G)`?

    # Let dl(H) denote the derived length of a group H.
    # Step 1: Find an upper bound for dl(hat(G)).
    # The condition that G_i / G_{i+1} is abelian implies [G_i, G_i] is a subgroup of G_{i+1}.
    # Since G = G_1, we have G' = [G, G] subset G_2.
    # By induction, G^(k) (the k-th term of the derived series) is a subgroup of G_{k+1}.
    # For k = n, we get G^(n) subset G_{n+1} = {1}. Thus, G^(n) is the trivial group.
    # This means dl(G) <= n.
    # The p-localization functor relates the derived series of G and `hat(G)` by
    # (hat(G))^(k) = hat(G^(k)).
    # If G^(k) = {1}, then hat(G^(k)) = {1}. So, dl(hat(G)) <= dl(G).
    # Combining these, we get dl(hat(G)) <= n.

    # Step 2: Show that the upper bound n is achievable.
    # The derived length of `hat(G)` is the smallest integer d such that G^(d) is a p'-group
    # (a group where every element's order is finite and not divisible by p).
    # To show the maximum derived length is n, we must find an example group G
    # satisfying the conditions for which G^(n-1) is not a p'-group.
    #
    # Let W_{n-1} be the (n-1)-fold iterated wreath product of the integers, Z.
    # W_{n-1} is a poly-infinite-cyclic group of derived length n-1. It is torsion-free,
    # so all factors in its characteristic series are p'-torsion-free.
    # Let Z_p be the cyclic group of order p.
    # Consider the group G = Z_p wr W_{n-1} (the wreath product).
    #
    # This group G fits the problem's criteria. It has a normal subgroup K
    # (the base group, a direct sum of copies of Z_p) such that G/K is isomorphic to W_{n-1}.
    # We can form a subnormal series for G of length n by taking the preimage of the
    # series for W_{n-1} and appending K at the end. The factors are p'-torsion-free
    # until the last one, which is K (an abelian group).
    #
    # The derived length of a wreath product A wr B is dl(A) + dl(B).
    # So, dl(G) = dl(Z_p) + dl(W_{n-1}) = 1 + (n-1) = n.
    # For a wreath product A wr B, the dl(B)-th derived subgroup is contained in the base group.
    # Here, G^(n-1) is a subgroup of the base group K.
    # Since dl(G) = n, G^(n-1) is not trivial. As a subgroup of K, it is a non-trivial p-group.
    # A non-trivial p-group is not a p'-group.
    #
    # Since G^(n-1) is not a p'-group, dl(hat(G)) must be at least n.
    # Combined with the upper bound dl(hat(G)) <= n, this proves that for this
    # specific G, dl(hat(G)) = n.
    # Therefore, the maximum possible derived length is n.
    answer_b_expression = "n"

    # The instruction "you still need to output each number in the final equation!"
    # is interpreted as outputting the symbolic expression 'n' for part (b).
    final_answer = f"(a) {answer_a}; (b) {answer_b_expression}"
    print(final_answer)

solve_group_theory_problem()