def absorption_cross_section_equations():
    """
    This function generates and prints the equations for the absorption
    cross-section for a chain of molecules under two different conditions.
    """

    # --- Explanation of Symbols ---
    explanation = """
    Here are the explanations for the symbols used in the equations below:

    σ(ω_L) : The absorption cross-section as a function of the laser's central frequency.
    ω_L    : The central angular frequency of the Gaussian laser pulse.
    ∝      : Proportional to. The equations omit fundamental constants for clarity.
    Σ_{i,a}  : A sum over all relevant transitions, from an initial occupied
             molecular state |i> to a final unoccupied molecular state |a>.
    N      : The total number of molecules in the chain.
    ω_{ia}   : The transition frequency of a single, isolated molecule for the
             i → a transition, given by (E_a - E_i)/ħ.
    μ_{ia}   : The transition dipole moment for the i → a transition in a
             single molecule. |μ_{ia}|² is the transition strength.
    τ      : A parameter related to the duration of the ultrashort laser pulse.
             A shorter pulse (smaller τ) has a broader frequency spectrum.
    J_{ia}   : The near-neighbor coupling energy (interaction strength) for an
             excitation of type i → a. It describes how an exciton on one
             molecule interacts with its neighbor.
    ħ      : The reduced Planck constant.

    The term 'exp(...)' represents the Gaussian lineshape function arising from the
    spectral profile of the ultrashort laser pulse.
    """
    print(explanation)

    # --- Case a) No interaction between molecules ---
    print("---------------------------------------------------------------------")
    print("a) Equation for non-interacting molecules:")
    print("---------------------------------------------------------------------")
    # In this case, the total absorption is the sum of N identical, independent absorbers.
    # The sum Σ_{n=1 to N} becomes a simple factor of N.
    equation_a = "σ(ω_L) ∝ Σ_{i,a} [ N ⋅ ω_{ia} ⋅ |μ_{ia}|² ⋅ exp(-(ω_{ia} - ω_L)² ⋅ τ²) ]"
    print(equation_a)
    print("\n   Each of the N molecules absorbs independently. The spectrum shows peaks\n"
          "   at each molecular transition frequency ω_{ia}, broadened by the laser pulse.")


    # --- Case b) Near-neighbor interaction considered ---
    print("\n---------------------------------------------------------------------")
    print("b) Equation with near-neighbor interactions (Frenkel excitons):")
    print("---------------------------------------------------------------------")
    # In this case, interactions create delocalized exciton states. A selection
    # rule dictates that only the K=0 state is optically active. This state has its
    # energy shifted by 2J and its transition strength enhanced by a factor of N.
    # The term for the transition dipole moment becomes |sqrt(N)μ_{ia}|² = N|μ_{ia}|².
    # The transition frequency is shifted to ω_{ia} + 2J_{ia}/ħ.
    equation_b = "σ(ω_L) ∝ Σ_{i,a} [ N ⋅ (ω_{ia} + 2J_{ia}/ħ) ⋅ |μ_{ia}|² ⋅ exp(-((ω_{ia} + 2J_{ia}/ħ) - ω_L)² ⋅ τ²) ]"
    print(equation_b)
    print("\n   Interactions create collective exciton states. The absorption is dominated\n"
          "   by a single, strong transition (for each i→a pair) which is shifted in\n"
          "   energy by the coupling J. This is known as a Frenkel exciton.")

# Execute the function to print the equations
absorption_cross_section_equations()
