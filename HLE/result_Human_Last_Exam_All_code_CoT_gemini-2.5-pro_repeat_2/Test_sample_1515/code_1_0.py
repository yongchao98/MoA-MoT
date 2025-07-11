def display_nsvz_condition_and_formula():
    """
    Explains the condition for the NSVZ beta function and displays its formula.
    """
    print("The NSVZ beta function is an exact result for the running of the gauge coupling in N=1 supersymmetric Yang-Mills (SYM) theories.")
    print("The non-renormalization theorems in SUSY are powerful statements about what quantum corrections are forbidden.")
    print("\nThe key connection between them lies in the concept of holomorphy. The gauge kinetic function is a holomorphic function of the complexified gauge coupling.")
    print("Non-renormalization theorems protect this structure. For the simple, exact form of the NSVZ beta function to hold, the regularization scheme used in calculations must respect this property.")
    print("\nTherefore, the exact condition is:")
    print("B. Regularization preserves holomorphy properties.")
    print("\nThis means that in a 'holomorphic' or 'NSVZ' scheme, the beta function is exactly related to the anomalous dimensions of the matter fields.")
    print("The formula is:")

    # Define the symbolic parts of the equation
    beta_g = "β(g)"
    g = "g"
    T_G = "T(G)"
    sum_T_r_gamma = "Σ_i T(r_i) * (1 - γ_i)"
    
    # Constants
    num_3 = "3"
    num_16 = "16"
    num_8 = "8"
    pi_sq = "π^2"

    # Construct the formula string by string
    numerator = f"{num_3}*{T_G} - {sum_T_r_gamma}"
    denominator = f"1 - ({g}^2 * {T_G}) / ({num_8}*{pi_sq})"
    
    # Print the final equation with all its parts
    print(f"\n{beta_g} = - ( {g}^3 / ({num_16}*{pi_sq}) ) * [ ({numerator}) / ({denominator}) ]\n")

    print("Where:")
    print(f"  {beta_g}: The beta function describing the running of the gauge coupling g.")
    print(f"  {g}: The gauge coupling constant.")
    print(f"  {T_G}: The Dynkin index for the adjoint representation of the gauge group.")
    print(f"  T(r_i): The Dynkin index for the i-th matter field representation.")
    print(f"  γ_i: The anomalous dimension of the i-th matter superfield.")

# Run the function to display the explanation and formula
display_nsvz_condition_and_formula()
<<<B>>>