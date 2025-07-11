import sympy

def derive_simplified_model():
    """
    This function derives and prints the simplified expression for synaptic efficacy dynamics.
    It uses the sympy library for symbolic mathematics to represent the variables and the equation.
    """
    # Define the symbols for the variables and parameters
    tau_w, w_i, u_i, v_i, rho = sympy.symbols('tau_w w_i u_i v_i rho')
    alpha, beta, eta = sympy.symbols('alpha beta eta')
    
    # Define the left-hand side of the equation
    lhs = tau_w * sympy.Symbol('dw_i/dt')
    
    # Define the right-hand side of the equation based on the derivation
    # Starting from: tau_w * dw_i/dt = alpha*Y + (beta-alpha)*B_i
    # With B_i = Y * (eta + M_i) / (1 + M_i)
    # This becomes: Y * (alpha + (beta-alpha)*(eta + M_i)/(1 + M_i))
    # Factoring out (beta-alpha): Y * (beta-alpha) * [alpha/(beta-alpha) + (eta+M_i)/(1+M_i)]
    # Substituting variables and rho: u_i * (beta-alpha) * [-rho + (eta+v_i)/(1+v_i)]
    
    # Create the expression for the right-hand side
    rhs = (beta - alpha) * u_i * (sympy.Symbol(f'(v_i + {eta})/(v_i + 1)') - rho)
    
    # Format the equation for printing
    # We will manually construct the string for a cleaner look
    final_equation = f"tau_w * dw_i/dt = (beta - alpha) * u_i * ( (v_i + eta) / (v_i + 1) - rho )"
    
    # Print the definitions of the simplified variables
    print("The simplified variables are defined as:")
    print("w_i: Synaptic efficacy W_i")
    print("u_i: Postsynaptic accumulator Y")
    print("v_i: Presynaptic accumulator M_i")
    print("")
    print("The constant rho is defined as:")
    print("rho = -alpha / (beta - alpha)")
    print("")
    
    # Print the final derived equation
    print("The derived expression for the dynamics of synaptic efficacy is:")
    print(final_equation)

if __name__ == '__main__':
    derive_simplified_model()
