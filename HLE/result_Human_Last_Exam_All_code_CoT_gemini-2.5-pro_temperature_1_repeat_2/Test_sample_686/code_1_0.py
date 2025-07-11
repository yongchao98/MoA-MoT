import sympy as sp

def display_susceptibility_equation():
    """
    This function displays the derived formula for the magnetic susceptibility
    of an Ising model on a sparse random graph.
    """
    # Define symbols for the equation using sympy for clear mathematical representation
    chi, beta, c, J, m_0, m_c, K, N = sp.symbols('chi beta c J m_0 m_c K N')
    tanh = sp.tanh

    # Define the expressions for N and K
    N_expr = beta * c * (1 - m_0**2) / (c - 1)
    K_expr = (tanh(beta * J) * (1 - m_c**2)) / (1 - (tanh(beta * J) * m_c)**2)
    
    # Define the final expression for chi
    chi_expr = N * (c - 1) * K / (1 - (c - 1) * K)

    print("The final expression for the magnetic susceptibility chi is:")
    print("==========================================================")
    # The final equation for chi
    print(f"chi = {chi_expr}")
    
    print("\nWhere the terms in the equation are defined as:")
    print("-----------------------------------------------")
    
    # The constant N
    print("The constant N:")
    print(f"N = {N_expr}")
    
    # The propagation factor K
    print("\nThe perturbation propagation factor K:")
    print(f"K = {K_expr}")
    
    print("\nNote:")
    print(" - 'beta' is the inverse temperature (1/k_B T).")
    print(" - 'c' is the connectivity of the graph.")
    print(" - 'J' is the homogeneous coupling constant.")
    print(" - 'm_0' is the site magnetization.")
    print(" - 'm_c' is the cavity magnetization.")
    print("The values of m_0 and m_c are determined by the model's self-consistency equations for a given external field B.")

# Execute the function to print the result
display_susceptibility_equation()

# For the final answer, we format the main equation as requested
# The following part will extract the components of the main formula for chi
# and format it as a string to be embedded in the final answer marker.
final_answer_str = f"chi = (N * (c - 1) * K) / (1 - (c - 1) * K)"
final_answer_str_simplified_numerator = f"beta * c * (1 - m_0**2) * K"
final_answer_str_simplified_denominator = f"1 - (c - 1) * K"

# We present the slightly more expanded version as the final answer,
# as it's more self-contained before introducing N.
final_answer = f"chi = ({final_answer_str_simplified_numerator}) / ({final_answer_str_simplified_denominator})"

print(f"\n<<<chi = (beta * c * (1 - m_0**2) * K) / (1 - (c - 1) * K)>>>")