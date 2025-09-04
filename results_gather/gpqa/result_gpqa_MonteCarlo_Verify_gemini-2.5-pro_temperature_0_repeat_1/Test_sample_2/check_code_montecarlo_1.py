import numpy as np
import math

def check_quantum_expectation_value():
    """
    This function checks the correctness of the given answer to the quantum mechanics problem.
    It calculates the expectation value of the operator 10*sigma_z + 5*sigma_x for the
    given state and compares it to the provided answer.
    """
    try:
        # The given answer is C, which corresponds to the value -0.7
        provided_answer_value = -0.7

        # Define the state vector |ψ⟩ = 0.5|↑⟩ + sqrt(3)/2|↓⟩
        # In the standard z-basis, |↑⟩ is represented by [1, 0] and |↓⟩ by [0, 1].
        c1 = 0.5
        c2 = math.sqrt(3) / 2
        psi = np.array([c1, c2])

        # First, check if the state vector is normalized.
        # The sum of the squares of the amplitudes must be 1.
        norm_squared = np.sum(np.abs(psi)**2)
        if not np.isclose(norm_squared, 1.0):
            return f"Incorrect: The state vector |ψ⟩ is not normalized. The sum of the squares of its coefficients is {norm_squared:.4f}, but it should be 1."

        # Define the Pauli matrices in the z-basis
        sigma_z = np.array([[1, 0], [0, -1]])
        sigma_x = np.array([[0, 1], [1, 0]])

        # Define the operator O = 10*sigma_z + 5*sigma_x
        operator_O = 10 * sigma_z + 5 * sigma_x

        # Calculate the expectation value ⟨ψ|O|ψ⟩.
        # This is calculated as ψ† * O * ψ, where ψ† is the conjugate transpose of ψ.
        # For a real vector psi, this is equivalent to psi.T @ O @ psi.
        # A more general way in numpy is np.vdot(psi, O @ psi).
        
        # Step 1: Apply the operator to the state -> O|ψ⟩
        O_psi = operator_O @ psi
        
        # Step 2: Take the inner product with the bra ⟨ψ| -> ⟨ψ|O|ψ⟩
        expectation_value = np.vdot(psi, O_psi)

        # The question asks for the value up to one decimal place.
        calculated_value_rounded = round(expectation_value.real, 1)

        # Compare the calculated value with the provided answer's value.
        if np.isclose(calculated_value_rounded, provided_answer_value):
            return "Correct"
        else:
            return (f"Incorrect: The calculated expectation value does not match the answer.\n"
                    f"The expectation value is ⟨ψ|O|ψ⟩ = 10⟨ψ|σ_z|ψ⟩ + 5⟨ψ|σ_x|ψ⟩.\n"
                    f"⟨ψ|σ_z|ψ⟩ = |c1|^2 - |c2|^2 = {c1**2 - c2**2} = -0.5.\n"
                    f"⟨ψ|σ_x|ψ⟩ = 2*Re(c1*c2) = 2 * {c1} * {c2:.4f} = {math.sqrt(3)/2:.4f} ≈ 0.866.\n"
                    f"Total expectation value = 10*(-0.5) + 5*({math.sqrt(3)/2:.4f}) = -5 + {2.5*math.sqrt(3):.4f} ≈ -0.67.\n"
                    f"Rounding to one decimal place gives {calculated_value_rounded}.\n"
                    f"The provided answer corresponds to {provided_answer_value}, which is consistent with the calculation.")

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Run the check
result = check_quantum_expectation_value()
print(result)