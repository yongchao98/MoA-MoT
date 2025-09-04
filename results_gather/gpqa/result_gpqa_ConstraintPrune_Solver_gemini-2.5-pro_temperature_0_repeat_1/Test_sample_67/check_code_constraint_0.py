import sympy as sp

def check_spin_oscillation_frequency():
    """
    This function symbolically calculates the oscillation frequency of the z-component
    of a spin-1/2 particle's magnetic moment and verifies it against the given answer.

    The problem describes a spin initially aligned with +Z, which then experiences
    a magnetic field in the +Y direction.
    """
    try:
        # Define the symbolic variables for the physical quantities.
        # gamma: gyromagnetic ratio
        # B: magnetic field strength
        # t: time
        # hbar: reduced Planck constant
        gamma, B, t, hbar = sp.symbols('gamma B t hbar', real=True, positive=True)
        i = sp.I  # Imaginary unit

        # --- Step 1: Define the initial state ---
        # The particle is initially spin-up along the Z-axis.
        # In the standard basis, |+z> is represented by the vector [1, 0].
        psi_0 = sp.Matrix([1, 0])

        # --- Step 2: Define the Hamiltonian ---
        # The magnetic field B is in the +Y direction. The Hamiltonian H = -μ · B.
        # Since μ = gamma * S, H = -gamma * S · (B * y_hat) = -gamma * B * Sy.
        # The spin-Y operator Sy is (hbar/2) * σy, where σy is the Pauli-Y matrix.
        sigma_y = sp.Matrix([[0, -i], [i, 0]])
        Sy = (hbar / 2) * sigma_y
        H = -gamma * B * Sy

        # --- Step 3: Calculate the time-evolved state ---
        # The state evolves according to the Schrödinger equation: |ψ(t)> = U(t)|ψ(0)>,
        # where the time-evolution operator U(t) = exp(-i*H*t/hbar).
        # We use sympy's matrix exponential function to compute U(t).
        U_t = sp.exp(-i * H * t / hbar)
        psi_t = sp.simplify(U_t * psi_0)

        # --- Step 4: Calculate the expectation value of the z-component of the moment ---
        # The observable is the z-component of the magnetic moment, μz = gamma * Sz.
        # We need to compute its expectation value: <μz(t)> = <ψ(t)| μz |ψ(t)>.
        # This is equal to gamma * <ψ(t)| Sz |ψ(t)>.
        # The spin-Z operator Sz is (hbar/2) * σz, where σz is the Pauli-Z matrix.
        sigma_z = sp.Matrix([[1, 0], [0, -1]])
        Sz = (hbar / 2) * sigma_z
        
        # Calculate <Sz(t)> = <ψ(t)| Sz |ψ(t)>. The .H attribute gives the conjugate transpose.
        exp_Sz_t_matrix = psi_t.H * Sz * psi_t
        # The result is a 1x1 matrix; we extract the scalar value.
        exp_Sz_t = sp.simplify(exp_Sz_t_matrix[0, 0])

        # The full expectation value of the z-moment is:
        exp_muz_t = gamma * exp_Sz_t
        
        # The derived expression is: gamma * (hbar/2) * cos(gamma*B*t)
        # This expression oscillates with time.

        # --- Step 5: Extract the oscillation frequency (ω) ---
        # The expression is of the form A*cos(ωt). The angular frequency ω is the
        # coefficient of 't' inside the cosine function's argument.
        
        # Find the cosine term in the expression
        cos_term = [arg for arg in exp_muz_t.args if isinstance(arg, sp.cos)]
        if not cos_term:
            return f"The expectation value <μz(t)> = {exp_muz_t} does not oscillate as a simple cosine function, which contradicts the problem's premise."
        
        # Extract the argument of the cosine function
        cos_argument = cos_term[0].args[0]
        
        # The frequency is the coefficient of time 't' in this argument.
        calculated_frequency = cos_argument.coeff(t)

        # --- Step 6: Verify against the provided answer ---
        # The provided answer is A, which corresponds to a frequency of gamma*B.
        answer_A_frequency = gamma * B
        
        if calculated_frequency == answer_A_frequency:
            return "Correct"
        else:
            return (f"Incorrect. The calculated oscillation frequency is {calculated_frequency}, "
                    f"but the frequency given in answer A is {answer_A_frequency}.")

    except Exception as e:
        return f"An error occurred during the verification process: {e}"

# Run the check
result = check_spin_oscillation_frequency()
print(result)