import textwrap

def explain_susy_qm_spectral_difference():
    """
    Explains the relationship between the spectra of two supersymmetric partner Hamiltonians
    and determines the maximum number of levels by which they can differ.
    """

    # Print the problem statement in a structured way
    print("Problem Analysis")
    print("================")
    print("We are given two Hamiltonians, H_0 and H_1, related by supersymmetry:")
    print("H_0 = L^+L - alpha")
    print("H_1 = L*L^+ - alpha")
    print("where L = d/dx - W(x) and L^+ = -d/dx - W(x).\n")

    # Step 1: Show the relationship between eigenstates
    print("Spectral Relationship")
    print("---------------------")
    explanation = """
    Let's see how the eigenstates of H_0 and H_1 are related.
    
    1. Mapping from H_0 to H_1:
    Let psi_0 be an eigenstate of H_0 with eigenvalue E.
    H_0 * psi_0 = E * psi_0
    (L^+L - alpha) * psi_0 = E * psi_0  =>  L^+L * psi_0 = (E + alpha) * psi_0
    
    Now, let's apply the operator L to this equation and see what happens to the state (L * psi_0):
    L * (L^+L * psi_0) = L * (E + alpha) * psi_0
    (L * L^+) * (L * psi_0) = (E + alpha) * (L * psi_0)
    
    Let's define a new state, psi_1 = L * psi_0. The equation becomes:
    (L * L^+) * psi_1 = (E + alpha) * psi_1
    
    Since H_1 = L*L^+ - alpha, we can write L*L^+ = H_1 + alpha. Substituting this in:
    (H_1 + alpha) * psi_1 = (E + alpha) * psi_1
    H_1 * psi_1 = E * psi_1
    
    This shows that if psi_0 is an eigenstate of H_0 with energy E, then psi_1 = L * psi_0 is an eigenstate of H_1 with the *same energy E*. This holds true as long as psi_1 is not zero.
    """
    print(textwrap.dedent(explanation))

    # Step 2: Analyze the special "unpaired" state
    print("The Special Case (Unpaired State)")
    print("-----------------------------------")
    explanation = """
    The mapping fails if L * psi_0 = 0. What is the energy of such a state?
    We go back to the H_0 eigenvalue equation:
    (L^+L - alpha) * psi_0 = E * psi_0
    
    If L * psi_0 = 0, this simplifies to:
    L^+(0) - alpha * psi_0 = E * psi_0
    -alpha * psi_0 = E * psi_0
    So, E = -alpha.
    
    This means that if there exists a normalizable state psi_0 such that L * psi_0 = 0, this state is an eigenstate of H_0 with energy -alpha. This state does not have a partner in the spectrum of H_1 because the mapping procedure gives zero.
    
    A symmetric argument holds for H_1. An eigenstate psi_1 of H_1 annihilated by L^+ (i.e., L^+ * psi_1 = 0) would be an unpaired eigenstate of H_1 with energy E = -alpha.
    """
    print(textwrap.dedent(explanation))

    # Step 3: Non-degeneracy and Conclusion
    print("Conclusion")
    print("----------")
    explanation = """
    A fundamental result from Sturm-Liouville theory states that for a one-dimensional Schrödinger operator on a finite interval (like H_0 and H_1 on [0,1]), the energy eigenvalues are non-degenerate.
    
    This implies:
    1. There can be at most one state (up to normalization) that satisfies L * psi_0 = 0.
    2. There can be at most one state that satisfies L^+ * psi_1 = 0.
    
    Furthermore, in supersymmetric quantum mechanics, it can be shown that at most one of these two potential unpaired states can be normalizable. It is not possible for both H_0 and H_1 to have an unpaired state.
    
    Therefore, the spectra of H_0 and H_1 are almost identical. All energy levels are paired and have the same value, with the possible exception of a single, non-degenerate level at energy E = -alpha that may exist in one spectrum but not the other.
    """
    print(textwrap.dedent(explanation))
    
    # Final Answer
    print("\nFinal Result")
    print("============")
    max_difference = 1
    print(f"The maximum number of levels of the spectrum of the Hamiltonians H_0 and H_1 that can differ is given by the equation:")
    print(f"maximum_difference = {max_difference}")


if __name__ == "__main__":
    explain_susy_qm_spectral_difference()
