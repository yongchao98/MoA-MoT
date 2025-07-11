import math

def display_nsvz_info():
    """
    Explains the condition for the NSVZ beta function and prints its formula.
    """
    print("The Novikov-Shifman-Vainshtein-Zakharov (NSVZ) beta function is an exact formula for the running of the gauge coupling 'g' in N=1 supersymmetric Yang-Mills theories.")
    print("It connects the evolution of the gauge coupling to the anomalous dimensions of the matter fields.")
    print("\n--- The Exact Condition ---")
    print("The validity of this exact, all-orders formula hinges on the choice of the renormalization scheme.")
    print("Specifically, the NSVZ formula holds in a scheme where the regularization procedure preserves the holomorphic structure of the theory.")
    print("Supersymmetry endows the theory with special properties related to holomorphy (quantities depending on complex variables in a complex-differentiable way). If the regularization scheme (e.g., a naive dimensional regularization) breaks this structure, the exact relation is violated.")
    print("Thus, the crucial condition is that the 'Regularization preserves holomorphy properties'.\n")

    print("--- NSVZ Beta Function Equation ---")
    # The user requested to output each number in the final equation.
    # We will print the formula as a formatted string.
    # beta(g) = - (g^3 / (16*pi^2)) * [3*C(G) - sum_i(T(R_i)*(1 - gamma_i))] / [1 - (g^2 / (8*pi^2))*C(G)]
    
    # We build the string representation to be clear.
    g_cubed = "g^3"
    sixteen_pi_squared = "16 * pi^2"
    
    three_C_G = "3 * C(G)"
    sum_term = "sum_i[T(R_i) * (1 - gamma_i)]"
    numerator = f"({three_C_G} - {sum_term})"
    
    one = "1"
    g_squared = "g^2"
    eight_pi_squared = "8 * pi^2"
    C_G = "C(G)"
    denominator = f"({one} - ({g_squared} / {eight_pi_squared}) * {C_G})"

    print("The equation is:")
    print(f"beta(g) = -({g_cubed} / {sixteen_pi_squared}) * [{numerator} / {denominator}]")
    
    print("\nWhere:")
    print("  beta(g): The beta function for the gauge coupling g.")
    print("  g:       The gauge coupling constant.")
    print("  C(G):    The quadratic Casimir of the adjoint representation.")
    print("  T(R_i):  The Dynkin index of the matter representation R_i.")
    print("  gamma_i: The anomalous dimension of the matter superfield i.")
    print("---------------------------------")

# Execute the function to print the information.
display_nsvz_info()