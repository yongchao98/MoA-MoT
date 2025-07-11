def solve():
    """
    This function formalizes the reasoning for the answers.

    (a) Is it true that sigma(a_j) = c_j a_{j-1}^* for a vertex j on the axis?
    Yes. The map sigma is expected to follow the geometry of the reflection g.
    g maps the source of a_j (e_j) to itself and the target (e_{j+1}) to e_{j-1}.
    So, sigma(a_j) must be an arrow from e_j to e_{j-1}.
    The fundamental arrow for this path is a_{j-1}^*. Thus, sigma(a_j) must be a scalar multiple of a_{j-1}^*.

    (b) Does sigma(a_j^*) = c_j^* a_j imply c_j^* = -mu_j^{-1} c_j?
    No. Following the same logic as in (a), g maps the source of a_j^* (e_{j+1}) to e_{j-1} and the target (e_j) to itself.
    So, sigma(a_j^*) must be an arrow from e_{j-1} to e_j, which is a multiple of a_{j-1}.
    The premise `sigma(a_j^*) = c_j^* a_j` can thus only hold if both sides are zero because a_j and a_{j-1} are linearly independent.
    This means the premise implies c_j^* = 0.
    The conclusion is `c_j^* = -mu_j^{-1} c_j`, which becomes `0 = -mu_j^{-1} c_j`, which is equivalent to `c_j = 0`.
    However, there is no reason why sigma(a_j^*) = 0 would imply sigma(a_j) = 0 (i.e., c_j = 0).
    So the implication is not necessary.

    (c) If sigma(a_i) is non-zero for an edge not on the axis, must lambda^2 mu_i mu_i^* = 1?
    No. This question likely refers to deformed preprojective algebras. For sigma to be an automorphism
    (e.g., of the form sigma(x) = lambda * g(x)), we can deduce two conditions if the deformation parameters lambda_i are not all zero:
    1. lambda^2 = 1
    2. (mu_k * mu_k^*)^2 = 1 for all k. Let C = mu_k * mu_k^*.
    The question is whether lambda^2 * C = 1.
    From our conditions, this is 1 * C = 1, so C = 1.
    But C^2 = 1 allows for C = -1.
    If C = -1, then lambda^2 * C = 1 * (-1) = -1, which is not 1.
    So the condition is not always true.
    """
    answer_a = "Yes"
    answer_b = "No"
    answer_c = "No"

    print(f"(a) [{answer_a}]; (b) [{answer_b.lower()}]; (c) [{answer_c.lower()}].")

solve()