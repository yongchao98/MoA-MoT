def check_nsvz_condition(regularization_scheme):
    """
    Checks if a given regularization scheme satisfies the condition for the NSVZ beta function to hold.
    The primary condition is the preservation of holomorphy.
    """
    print(f"Analyzing the scheme: '{regularization_scheme['name']}'...")

    # The exact condition for the NSVZ beta function to match non-renormalization theorems
    # is that the regularization scheme must preserve the holomorphy properties of the theory.
    if regularization_scheme['preserves_holomorphy']:
        print("✅ Condition met: The regularization scheme preserves holomorphy.")
        print("This allows the one-loop exact beta function for the holomorphic coupling to be related directly to the physical beta function.")
        print("As a result, the NSVZ beta function formula is expected to hold.")

        # The NSVZ beta function is given by:
        # β(g) = - [g^3 / (16 * π^2)] * [3*T(adj) - Σ T(R_i)*(1 - γ_i)] / [1 - g^2*T(adj)/(8*π^2)]
        # We will print the numerical constants and exponents from this famous equation.
        print("\nSymbolic representation of the NSVZ beta function equation:")
        print("β(g) ∝ -g^3 * (3*T_adj - ΣT_R(1-γ)) / (1 - C*g^2)")
        print("\nKey numerical values and exponents in the equation:")
        print(f"Power of the gauge coupling 'g' in the numerator: 3")
        print(f"Coefficient of the adjoint Casimir T(adj) in the numerator: 3")
        print(f"Coefficient of the matter term (1-γ): 1")
        print(f"Power of the gauge coupling 'g' in the denominator term: 2")
        print(f"Numerical factor in the denominator of the main prefactor (16*π^2): 16")
        print(f"Numerical factor in the denominator of the g^2 term (8*π^2): 8")

    else:
        print("❌ Condition not met: The regularization scheme does not preserve holomorphy.")
        print("The standard NSVZ formula is not guaranteed to hold; scheme-dependent corrections are expected.")

# We model two types of schemes.
# Dimensional Reduction (DRED) is designed to preserve supersymmetry and thus holomorphy.
dred_scheme = {
    "name": "Dimensional Reduction",
    "preserves_holomorphy": True
}

# Let's execute the check with the appropriate scheme.
check_nsvz_condition(dred_scheme)
