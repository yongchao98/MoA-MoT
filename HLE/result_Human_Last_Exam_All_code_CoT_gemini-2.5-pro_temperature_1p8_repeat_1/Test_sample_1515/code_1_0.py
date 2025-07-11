import math

def print_nsvz_beta_function():
    """
    This function prints the formula for the NSVZ beta function.

    The NSVZ beta function provides an exact expression for the running of the gauge
    coupling in N=1 supersymmetric Yang-Mills theories coupled to matter.
    Its exactness relies on the choice of a regularization scheme that
    preserves the holomorphic structure of the theory.

    The formula is:
    β(g) = - [g^3 / (16*π^2)] * [3*C2(G) - T(R)*(1 - γ_m)] / [1 - (g^2 / (8*π^2))*C2(G)]
    """
    
    # Placeholder symbolic variables for the equation's components
    g_sq = "g^2"
    g_cubed = "g^3"
    pi_sq = "π^2"
    C2_G = "C2(G)"  # Dynkin index of the adjoint representation, T(adj)
    T_R = "T(R)"    # Dynkin index of the matter representation
    gamma_m = "γ_m"   # Anomalous dimension of the matter superfield

    # Building the numerator and denominator strings
    # Numerator: 3*C2(G) - T(R)*(1 - γ_m)
    numerator_expr = f"{3}*{C2_G} - {T_R}*(1 - {gamma_m})"
    
    # Denominator: 1 - (g^2 / (8*π^2))*C2(G)
    denominator_expr = f"1 - ({g_sq} / ({8}*{pi_sq}))*{C2_G}"

    # Construct the final equation string
    # We output each number in the equation explicitly.
    print("The NSVZ Beta Function Equation:")
    print(f"β(g) = - ( {g_cubed} / ({16}*{pi_sq}) ) * ( ( {numerator_expr} ) / ( {denominator_expr} ) )")
    print("\nWhere:")
    print("  g: Gauge coupling constant")
    print("  γ_m: Anomalous dimension of the matter chiral superfield")
    print("  C2(G): Second Casimir invariant for the adjoint representation")
    print("  T(R): Dynkin index for the matter representation R")

    print("\nThe condition for this formula to be exact to all loops is that the regularization scheme must preserve the holomorphy properties of the supersymmetric theory.")

# Execute the function to print the equation
print_nsvz_beta_function()
