import sys

def solve_exponent_order():
    """
    This script explains and determines the order in the coupling constant 'u'
    at which the critical exponent 'nu' receives its first non-vanishing contribution
    within the framework of the epsilon-expansion for phi^4 theory.
    """

    # 1. The relationship between nu and the anomalous dimension of phi^2.
    #    The anomalous dimension is a function of the coupling constant u, denoted gamma_2(u).
    print("Step 1: The foundational equation from Renormalization Group theory.")
    print("The critical exponent ν is related to the anomalous dimension of the φ² operator, γ₂(u):")
    # Using unicode for clarity in the output.
    if sys.stdout.encoding.lower().startswith('utf'):
        print("  1 / ν = 2 - γ₂(u)")
    else:
        print("  1 / nu = 2 - gamma_2(u)")
    print("-" * 50)

    # 2. Perturbative expansion of the anomalous dimension.
    #    Loop calculations show that the expansion of gamma_2(u) starts with a linear term.
    print("Step 2: Perturbative expansion of the anomalous dimension.")
    print("The function γ₂(u) is calculated as a power series in the coupling constant u.")
    print("The one-loop contribution is non-zero, making the expansion start at order u¹:")
    if sys.stdout.encoding.lower().startswith('utf'):
        print("  γ₂(u) = A·u¹ + B·u² + O(u³)")
    else:
        print("  gamma_2(u) = A*u^1 + B*u^2 + O(u^3)")
    print("  (where 'A' is a non-zero constant from the one-loop calculation)")
    print("-" * 50)


    # 3. Deriving the correction to nu.
    #    Substituting the series for gamma_2(u) into the equation for nu.
    print("Step 3: Calculating the correction for ν.")
    print("Substituting the series for γ₂(u) into the main equation:")
    if sys.stdout.encoding.lower().startswith('utf'):
        print("  1 / ν = 2 - (A·u¹ + ...)")
    else:
        print("  1 / nu = 2 - (A*u^1 + ...)")
    print("\nSolving for ν and expanding for small u gives:")
    if sys.stdout.encoding.lower().startswith('utf'):
        print("  ν = 1/2 + (A/4)·u¹ + O(u²)")
    else:
        print("  nu = 1/2 + (A/4)*u^1 + O(u^2)")
    print("-" * 50)


    # 4. Final Conclusion.
    order = 1
    print("Conclusion:")
    print("The first correction to the mean-field value of ν=1/2 is the term proportional to u¹.")
    print(f"Thus, the initial non-vanishing contribution appears at order {order} in the coupling constant u.")

    print("\nThe question asks for the specific order N in the relation: Correction(ν) ∝ u^N")
    print("The number in the final equation is:")
    print(f"N = {order}")

solve_exponent_order()
