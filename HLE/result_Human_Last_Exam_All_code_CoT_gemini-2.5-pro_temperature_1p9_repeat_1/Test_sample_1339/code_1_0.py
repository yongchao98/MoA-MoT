import math

def solve_group_theory_problem():
    """
    This function solves the given group theory problem.

    The reasoning is as follows:
    1.  The notation `G = G_1 \triangleleft G_2 \triangleleft \dots \triangleleft G_n \triangleleft G_{n+1} = \{1\}` is interpreted as a chain of subgroups where G_i is a normal subgroup of G_{i+1}.
    2.  `A \triangleleft B` implies `A` is a subset of `B`. Thus, `G \subseteq G_2 \subseteq \dots \subseteq G_n \subseteq G_{n+1}`.
    3.  Given `G_{n+1} = \{1\}`, it follows that `G` must be a subgroup of the trivial group, meaning `G = \{1\}`.
    4.  For G = \{1\}, any p-nonsingular system has a solution in G itself. Thus, the minimal overgroup `\hat{G}` is G itself, i.e., `\hat{G} = \{1\}`. This group is unique up to isomorphism.
    5.  The derived length of the trivial group `\hat{G}=\{1\}` is 0.
    """
    
    # Part (a)
    # Based on the reasoning, a unique minimal group `\hat{G}` exists.
    answer_a = "Yes"
    
    # Part (b)
    # The derived length of `\hat{G} = \{1\}` is 0. This is the maximum possible value since G must be trivial.
    derived_length = 0
    answer_b = derived_length

    # The problem asks for the answer in the form: (a) [Yes/No]; (b) [Expression].
    # The code below prints the final answer. The number in the expression is explicitly given.
    print(f"(a) {answer_a}; (b) {answer_b}")

solve_group_theory_problem()