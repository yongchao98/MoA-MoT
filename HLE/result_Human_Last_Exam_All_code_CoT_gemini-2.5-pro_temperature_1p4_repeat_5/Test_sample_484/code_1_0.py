def derive_simplified_model():
    """
    This script performs a steady-state analysis on the provided biophysical model
    to derive a simplified equation for synaptic efficacy dynamics.
    """
    
    print("Derivation of the simplified synaptic plasticity model\n")
    print("="*50)
    print("Step 1: Steady-State Analysis of Fast Variables")
    print("="*50)
    print("We assume the variables M, Y, P, and B change much faster than the synaptic weights W.")
    print("Thus, we can analyze their steady state by setting their time derivatives to zero.")
    print("We also replace the spike train input x_i(t) with its average firing rate, ν_i.\n")

    print("1.1) Presynaptic MMP9 (M_i):")
    print("Original equation: τ_M * dM_i/dt = -M_i + φ*x_i(t)")
    print("Setting dM_i/dt = 0 and x_i(t) = ν_i gives:")
    print("=> 0 = -M_i + φ*ν_i")
    print("=> M_i = φ*ν_i\n")

    print("1.2) Postsynaptic Calcium (Y):")
    print("Original equation: τ_Y * dY/dt = -Y + Σ_j(w_j * x_j(t))")
    print("Setting dY/dt = 0 and x_j(t) = ν_j gives:")
    print("=> 0 = -Y + Σ_j(w_j * ν_j)")
    print("=> Y = Σ_j(w_j * ν_j)\n")

    print("1.3) proBDNF (P_i) and BDNF (B_i):")
    print("Original equations:")
    print("  τ_P * dP_i/dt = -P_i + (1-η)*Y - M_i*P_i")
    print("  τ_P * dB_i/dt = -B_i + η*Y + M_i*P_i")
    print("Setting the derivatives to zero gives:")
    print("=> P_i * (1 + M_i) = (1-η)*Y  =>  P_i = (1-η)*Y / (1 + M_i)")
    print("=> B_i = η*Y + M_i*P_i")
    print("Substituting P_i into the expression for B_i:")
    print("=> B_i = η*Y + M_i*[(1-η)*Y / (1 + M_i)] = Y * [η + M_i*(1-η)/(1+M_i)]")
    print("=> B_i = Y * [η(1+M_i) + M_i(1-η)] / (1 + M_i)")
    print("=> B_i = Y * (η + M_i) / (1 + M_i)\n")
    
    print("="*50)
    print("Step 2: Substitute into the Synaptic Efficacy Equation")
    print("="*50)
    print("The dynamics of the synaptic efficacy w_i (same as W_i) are given by:")
    print("τ_W * dw_i/dt = α*P_i + β*B_i\n")
    print("Substitute the steady-state expressions for P_i and B_i:")
    print("=> τ_W * dw_i/dt = α*[(1-η)*Y / (1+M_i)] + β*[Y*(η+M_i) / (1+M_i)]")
    print("Factor out Y / (1+M_i):")
    print("=> τ_W * dw_i/dt = [Y / (1+M_i)] * [α*(1-η) + β*(η+M_i)]")
    print("=> τ_W * dw_i/dt = [Y / (1+M_i)] * [α - αη + βη + βM_i]\n")

    print("="*50)
    print("Step 3: Introduce Simplified Variables (u_i, v_i, ρ)")
    print("="*50)
    print("We define new variables to represent the key accumulators and parameters:")
    print("  - Presynaptic accumulator (v_i): Represents presynaptic activity.")
    print("    Let v_i = M_i = φ*ν_i")
    print("  - Postsynaptic accumulator (u_i): Represents postsynaptic activity.")
    print("    Let u_i = Y = Σ_j(w_j*ν_j)")
    print("  - LTD/LTP ratio constant (ρ): Represents the balance of depression and potentiation.")
    print("    Let ρ = -α/β, which means α = -ρ*β\n")
    
    print("="*50)
    print("Step 4: Final Derivation")
    print("="*50)
    print("Substitute α = -ρ*β into the equation from Step 2:")
    print("τ_W * dw_i/dt = [Y / (1+M_i)] * [(-ρβ) - (-ρβ)η + βη + βM_i]")
    print("Factor out β:")
    print("=> τ_W * dw_i/dt = [β*Y / (1+M_i)] * [-ρ + ρη + η + M_i]")
    print("Group terms:")
    print("=> τ_W * dw_i/dt = [β*Y / (1+M_i)] * [M_i - ρ + η(1+ρ)]\n")
    print("Finally, replace Y and M_i with u_i and v_i:")
    print("The final expression for τ_w * dw_i/dt is:\n")
    
    final_equation = "τ_W * dw_i/dt = (β * u_i / (1 + v_i)) * (v_i - ρ + η*(1 + ρ))"
    print(final_equation)

if __name__ == '__main__':
    derive_simplified_model()
    # The final derived expression is the answer to the user's request.
    # It is derived step-by-step above and printed as the final output.
    final_expression_string = "τ_W*dw_i/dt = (β*u_i/(1 + v_i))*(v_i - ρ + η*(1 + ρ))"
    # To conform to the output format, the final answer is provided below.
    # print(f"<<<{final_expression_string}>>>")