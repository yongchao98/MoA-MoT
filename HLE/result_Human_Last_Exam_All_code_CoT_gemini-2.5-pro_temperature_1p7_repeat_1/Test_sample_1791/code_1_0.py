import sympy

def solve():
    """
    This function symbolically derives the effective electron-electron interaction
    Hamiltonian mediated by phonons using the method of completing the square.
    """
    # 1. Define symbolic variables
    # g (coupling strength) and m (mass) are real, positive constants.
    g, m = sympy.symbols('g m', real=True, positive=True)
    
    # rho_q is the electron density operator. Its Hermitian conjugate rho_q^dagger is rho_{-q}.
    rho_q = sympy.Symbol('rho_q')
    rho_mq = sympy.Symbol('rho_-q')
    
    # w_q is the phonon frequency. It's assumed to be a positive real value.
    omega_q = sympy.Symbol('w_q', real=True, positive=True)
    
    # We consider a 3D system where the wavevector q = (q_x, q_y, q_z).
    q_components = sympy.symbols('q_x q_y q_z', real=True)
    
    # Define q^2 as a symbol for |q|^2 for a cleaner final expression.
    q_squared = sympy.Symbol('q^2')

    # 2. Derive the interaction coefficient U(q)
    # The contribution from a single polarization j with wavevector component q_j
    q_j = sympy.Symbol('q_j', real=True)
    
    # The coupling constant C_{q,j} given in the problem
    C_qj = g * sympy.I * q_j / sympy.sqrt(2 * m * omega_q)

    # The effective interaction H_eff = - sum_{q,j} |C_{q,j}|^2/w_q * rho_q * rho_{-q}
    # First, calculate the squared magnitude of C_{q,j}
    C_qj_abs_sq = (C_qj * sympy.conjugate(C_qj)).simplify()

    # Now, find the interaction term for a single mode (q,j)
    # This is the coefficient of rho_q * rho_{-q}
    U_qj = -C_qj_abs_sq / omega_q

    # 3. Sum over polarizations j (i.e., j = x, y, z)
    # The terms only differ by the component of q, so summing them
    # involves replacing q_j^2 with q_x^2 + q_y^2 + q_z^2.
    total_U_q_symbolic = sum(U_qj.subs(q_j, qc) for qc in q_components)
    
    # For a more compact and general notation, we replace the sum of squares
    # of the components with the symbol for the squared magnitude of the vector.
    total_U_q = total_U_q_symbolic.subs(sum(qc**2 for qc in q_components), q_squared)

    # 4. Print the final results
    print("The effective electron-electron interaction Hamiltonian for a given wavevector q has the form:")
    print("H_eff(q) = U(q) * rho_q * rho_-q")
    print("\nwhere rho_q is the electron density operator.")
    print("The effective interaction potential U(q) is derived as:")
    
    # Pretty print the final symbolic expression for U(q)
    sympy.pprint(total_U_q)

    print("\nThus, the final equation for the effective Hamiltonian at a specific q is:")
    
    # Construct and print the full equation H_eff(q) = ...
    H_eff_q_symbol = sympy.Symbol('H_eff(q)')
    final_equation = sympy.Eq(H_eff_q_symbol, total_U_q * rho_q * rho_mq)
    sympy.pprint(final_equation, use_unicode=True)

solve()