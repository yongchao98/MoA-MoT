import sympy
from sympy import symbols, sin, cos, sqrt, Rational, expand, simplify

def check_quantum_energy_spectrum():
    """
    This function verifies the derivation of the energy spectrum for the given potential.
    It follows these steps:
    1. Defines the potential in polar coordinates.
    2. Converts the potential to Cartesian coordinates.
    3. Identifies the system as two independent harmonic oscillators and finds their frequencies.
    4. Calculates the total energy spectrum by summing the individual energies.
    5. Compares the derived result with the expression from the provided answer.
    """
    try:
        # Define symbolic variables for the calculation
        r, theta, x, y, k, m, hbar = symbols('r theta x y k m hbar', real=True, positive=True)
        n_x, n_y = symbols('n_x n_y', integer=True, nonneg=True)

        # Step 1: Define the potential in polar coordinates
        # V(r, θ) = 1/2 kr^2 + 3/2 kr^2 cos^2(θ)
        V_polar = Rational(1, 2) * k * r**2 + Rational(3, 2) * k * r**2 * cos(theta)**2

        # Step 2: Convert the potential to Cartesian coordinates
        # Use the transformations: x = r*cos(θ) and r^2 = x^2 + y^2
        # We can substitute r^2*cos^2(θ) with x^2, and then r^2 with x^2 + y^2
        V_cartesian = V_polar.subs(r**2 * cos(theta)**2, x**2).subs(r**2, x**2 + y**2)
        V_cartesian_simplified = simplify(V_cartesian)

        # The expected potential is V(x,y) = 2*k*x^2 + 1/2*k*y^2
        expected_V = 2*k*x**2 + Rational(1, 2)*k*y**2
        if simplify(V_cartesian_simplified - expected_V) != 0:
            return f"Incorrect potential conversion. Expected {expected_V}, but got {V_cartesian_simplified}."

        # Step 3: Determine effective spring constants and angular frequencies
        # The potential V(x,y) = 1/2 * kx_eff * x^2 + 1/2 * ky_eff * y^2
        # Comparing with V(x,y) = 2*k*x^2 + 1/2*k*y^2, we find:
        # 1/2 * kx_eff = 2*k  => kx_eff = 4*k
        # 1/2 * ky_eff = 1/2*k => ky_eff = k
        kx_eff = 4 * k
        ky_eff = k

        # The angular frequency is ω = sqrt(k_eff / m)
        omega_x = sqrt(kx_eff / m)
        omega_y = sqrt(ky_eff / m)

        # Step 4: Calculate the total energy spectrum
        # The energy for a 1D QHO is E_n = (n + 1/2)ħω
        # The total energy is E = E_x + E_y
        E_total = (n_x + Rational(1, 2)) * hbar * omega_x + \
                  (n_y + Rational(1, 2)) * hbar * omega_y
        
        # Simplify the final expression by expanding and collecting terms
        E_final_derived = simplify(expand(E_total))

        # Step 5: Compare the derived result with the provided answer's expression
        # The provided answer is B, which corresponds to the expression:
        # E = (2*n_x + n_y + 3/2) * hbar * sqrt(k/m)
        answer_B_expr = (2*n_x + n_y + Rational(3, 2)) * hbar * sqrt(k/m)

        # Check if the derived energy matches the expression for option B
        if simplify(E_final_derived - answer_B_expr) == 0:
            return "Correct"
        else:
            return (f"The derived energy spectrum is {E_final_derived}, which does not match the "
                    f"expression for the given answer B: {answer_B_expr}. The final simplification step in the provided answer's reasoning is correct, but the code check reveals a discrepancy in the overall logic or the provided options.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_quantum_energy_spectrum()
print(result)