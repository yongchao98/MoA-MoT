def derive_simplified_model():
    """
    This function prints the step-by-step derivation of the simplified model
    for synaptic plasticity.
    """

    print("--- Biophysical Model Simplification ---\n")
    print("Goal: Derive an expression for tau_w * dw_i/dt in terms of u_i, v_i, and rho.")
    print("This is achieved by performing a steady-state analysis on the fast variables (M, Y, P, B).\n")

    print("Step 1: Define the simplified variables (as per the prompt).")
    print("Let w_i be the synaptic efficacy.")
    print("Let v_i be the presynaptic accumulator.")
    print("Let u_i be the postsynaptic accumulator.\n")

    print("Step 2: Find the steady-state (ss) of the fast variables.\n")
    
    print("   a) For presynaptic MMP9 (M_i):")
    print("      Original equation: tau_M * dM_i/dt = -M_i + phi * x_i")
    print("      Setting dM_i/dt = 0 gives: M_i_ss = phi * x_i")
    print("      We define the presynaptic accumulator v_i as this steady-state value:")
    print("      v_i = M_i_ss = phi * x_i\n")

    print("   b) For postsynaptic calcium (Y):")
    print("      Original equation: tau_Y * dY/dt = -Y + sum_j(w_j * x_j)")
    print("      Setting dY/dt = 0 gives: Y_ss = sum_j(w_j * x_j)")
    print("      We define the postsynaptic accumulator u_i as this shared steady-state value:")
    print("      u_i = Y_ss = sum_j(w_j * x_j)\n")

    print("   c) For proBDNF (P_i) and BDNF (B_i):")
    print("      Original equations:")
    print("      tau_P * dP_i/dt = -P_i + (1-eta)*Y - M_i*P_i")
    print("      tau_P * dB_i/dt = -B_i + eta*Y + M_i*P_i")
    print("      Setting derivatives to zero, we get a system of two algebraic equations:")
    print("      P_i * (1 + M_i) = (1-eta)*Y  =>  P_i_ss = (1-eta)*Y / (1 + M_i)")
    print("      B_i = eta*Y + M_i*P_i")
    print("      Substituting P_i_ss into the equation for B_i gives:")
    print("      B_i_ss = eta*Y + M_i*[(1-eta)*Y / (1 + M_i)]")
    print("      B_i_ss = Y * [eta + M_i*(1-eta)/(1+M_i)] = Y * [ (eta*(1+M_i) + M_i*(1-eta)) / (1+M_i) ]")
    print("      B_i_ss = Y * [ (eta + eta*M_i + M_i - eta*M_i) / (1+M_i) ] = Y * [ (eta + M_i) / (1+M_i) ]\n")

    print("Step 3: Substitute the steady-state expressions into the weight dynamics equation.\n")
    print("      The weight dynamics are governed by:")
    print("      tau_w * dw_i/dt = alpha*P_i + beta*B_i")
    print("      Substitute the steady-state expressions for P_i and B_i:")
    print("      tau_w * dw_i/dt = alpha*[ (1-eta)*Y / (1+M_i) ] + beta*[ Y*(eta+M_i) / (1+M_i) ]")
    print("      Factor out Y / (1+M_i):")
    print("      tau_w * dw_i/dt = (Y / (1+M_i)) * [ alpha*(1-eta) + beta*(eta+M_i) ]\n")

    print("Step 4: Introduce the simplified variables u_i and v_i and rearrange.\n")
    print("      Substitute Y with u_i and M_i with v_i:")
    print("      tau_w * dw_i/dt = (u_i / (1+v_i)) * [ alpha*(1-eta) + beta*eta + beta*v_i ]")
    print("      To simplify, we factor out beta from the term in the brackets:")
    print("      tau_w * dw_i/dt = (beta*u_i / (1+v_i)) * [ (alpha*(1-eta) + beta*eta)/beta + v_i ]")
    print("      Rearranging to get a (v_i - threshold) form:")
    print("      tau_w * dw_i/dt = (beta*u_i / (1+v_i)) * [ v_i - (- (alpha*(1-eta) + beta*eta)/beta) ]\n")
    
    print("Step 5: Define the constant rho and state the final expression.\n")
    print("      We define the constant rho as the threshold term:")
    print("      rho = - (alpha*(1-eta) + beta*eta) / beta")
    print("      This can also be written as: rho = - (alpha/beta)*(1-eta) - eta\n")
    print("      Substituting rho into the equation gives the final simplified expression:")
    print("      =======================================================")
    final_equation = "      tau_w * dw_i/dt = (beta * u_i * (v_i - rho)) / (1 + v_i)"
    print(final_equation)
    print("      =======================================================\n")
    
if __name__ == '__main__':
    derive_simplified_model()
