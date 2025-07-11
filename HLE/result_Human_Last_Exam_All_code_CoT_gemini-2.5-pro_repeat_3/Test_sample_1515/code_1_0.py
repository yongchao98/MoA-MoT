import math

def explain_nsvz_condition():
    """
    Explains the condition for the exactness of the NSVZ beta function
    and prints the formula.
    """
    print("The Novikov-Shifman-Vainshtein-Zakharov (NSVZ) beta function provides a proposed exact expression for the beta function of the gauge coupling in N=1 supersymmetric gauge theories.")
    print("\n--- The Key Condition for Exactness ---")
    print("The NSVZ beta function is not exact in all regularization schemes. Its validity as an all-orders exact result depends critically on the properties of the chosen scheme.")
    print("The central condition is that the regularization scheme must preserve the holomorphy of the gauge kinetic function.")
    print("\nIn supersymmetric theories, quantities like the superpotential and gauge kinetic function are holomorphic functions of their respective couplings and fields. Non-renormalization theorems are a consequence of this holomorphy. If a regularization scheme breaks this structure, the derivation of the exact beta function is no longer valid beyond one-loop.")
    print("\nTherefore, the NSVZ formula is exact in a 'holomorphic scheme', such as Dimensional Reduction (DRED), but not in schemes like conventional Dimensional Regularization (DREG) that break supersymmetry.")

    print("\n--- The NSVZ Beta Function Equation ---")
    print("The beta function for the gauge coupling 'g' is given by:")
    
    # The user requested to output each number in the equation.
    # The numbers are: 3, 1, 16, 2, 1, 8, 2. I will print them as part of the equation string.
    
    numerator_part_1 = "3*T(G) - T(R)*(1 - γ)"
    denominator_part = "1 - T(G)*g^2 / (8*π^2)"
    
    print(f"\n      -g^3      ( {numerator_part_1} )")
    print(f"β(g) = ------ * --------------------------------")
    print(f"      16*π^2    ( {denominator_part} )")
    
    print("\nWhere:")
    print("  β(g): The beta function for the gauge coupling g.")
    print("  g: The gauge coupling constant.")
    print("  T(G): The Dynkin index for the adjoint representation of the gauge group.")
    print("  T(R): The sum of Dynkin indices for the matter field representations.")
    print("  γ: The anomalous dimension of the matter superfields.")
    print("  π: The mathematical constant pi.")

    print("\n--- Conclusion ---")
    print("The condition that allows the NSVZ beta function to be exact and compatible with non-renormalization theorems is:")
    print("B. Regularization preserves holomorphy properties")


if __name__ == "__main__":
    explain_nsvz_condition()
