import sympy as sp

def derive_simplified_plasticity_rule():
    """
    Performs a steady-state analysis on the biophysical model to derive a simplified
    three-factor learning rule.
    """
    # Define symbolic variables and parameters
    N = sp.Symbol('N', integer=True, positive=True)
    tau_W, phi, eta, alpha, beta = sp.symbols('tau_W phi eta alpha beta', real=True)
    
    # Use IndexedBase for indexed variables like firing rates (nu) and weights (w)
    i, j = sp.symbols('i j', cls=sp.Idx)
    nu = sp.IndexedBase('nu', real=True, positive=True)
    w = sp.IndexedBase('w', real=True)

    # --- Step 1 & 2: Steady-state solutions for fast variables ---
    # In steady-state, we use the mean firing rate nu_i instead of the spike train x_i(t).
    # M_i is the steady-state value of MMP9 at synapse i.
    M_i_ss = phi * nu[i]

    # Y_ss is the steady-state value of the shared postsynaptic calcium.
    # The sum is over all synapses j from 1 to N.
    Y_ss = sp.Sum(w[j] * nu[j], (j, 1, N))

    # P_i_ss and B_i_ss are the steady-state values of proBDNF and BDNF.
    P_i_ss = (1 - eta) * Y_ss / (1 + M_i_ss)
    B_i_ss = (eta * Y_ss + M_i_ss * P_i_ss)
    
    # --- Step 3: Substitute into the weight dynamics equation ---
    # The expression for the change in weight w_i.
    tau_W_dwidt = alpha * P_i_ss + beta * B_i_ss
    
    # Simplify the expression
    tau_W_dwidt_simplified = sp.simplify(tau_W_dwidt)
    # The result is: Y_ss*(alpha*(1 - eta) + beta*(eta + M_i_ss))/(1 + M_i_ss)

    # --- Step 4: Rearrange and define new variables ---
    # We want to get the form u_i * (v_i - rho)
    
    # Define v_i as the presynaptic accumulator, which is the steady-state MMP9 level.
    v_i_def = M_i_ss
    v_i_sym = sp.Symbol('v_i')

    # Define rho as the modification threshold.
    # We factor out beta from the numerator to find the expression for rho.
    # Numerator = alpha*(1 - eta) + beta*eta + beta*M_i_ss
    #           = beta * (M_i_ss - (- (alpha*(1-eta) + beta*eta)/beta))
    rho_def = - (alpha*(1 - eta) + beta*eta) / beta
    rho_sym = sp.Symbol('rho')
    
    # Define u_i as the postsynaptic accumulator.
    u_i_def = (beta * Y_ss) / (1 + M_i_ss)
    u_i_sym = sp.Symbol('u_i')

    # The final expression in the simplified form.
    final_expr = u_i_sym * (v_i_sym - rho_sym)
    
    # --- Print the derivation and result ---
    print("The steady-state analysis reduces the system to a simpler form.")
    print("The resulting expression for the change in synaptic efficacy is derived as follows:")

    print("\n1. The equation for the dynamics of synaptic efficacy is:")
    print(f"   τ_W * dW_i/dt = α*P_i + β*B_i")

    print("\n2. After substituting the steady-state values for P_i and B_i, we get:")
    # Use pretty print for better formatting of the sum
    pretty_expr = sp.printing.pretty(tau_W_dwidt_simplified, use_unicode=False)
    print(f"   τ_W * dW_i/dt = {pretty_expr}")

    print("\n3. This expression can be rewritten in a compact three-factor form by defining:")
    print(f"   - Presynaptic accumulator: v_i = {v_i_def}")
    print(f"   - Postsynaptic accumulator: u_i = {u_i_def}")
    print(f"   - Constant threshold: ρ = {sp.simplify(rho_def)}")

    print("\n4. The final simplified expression for the synaptic weight dynamics is:")
    print(f"   τ_W * dW_i/dt = u_i * (v_i - ρ)")
    
    # The final answer in the requested format.
    # The question asks for the expression for tau_w*dw_i
    final_answer_str = f"{u_i_sym}*(v_i - {rho_sym})"
    print(f"\n<<<{final_expr}>>>")


if __name__ == '__main__':
    derive_simplified_plasticity_rule()