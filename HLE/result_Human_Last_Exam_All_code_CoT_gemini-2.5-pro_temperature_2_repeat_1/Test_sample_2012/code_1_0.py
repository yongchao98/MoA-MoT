import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum mechanics problem based on the rules of Universe U.
    """
    # Step 1: Define the observable operator
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 2: Calculate eigenvalues and eigenvectors.
    # np.linalg.eigh sorts eigenvalues in ascending order.
    eigenvalues, eigenvectors = np.linalg.eigh(O)

    # For a 2x2 matrix, the eigenvalues will be sorted smallest to largest.
    # Therefore, eigenvalues[0] is the second-largest (smaller) eigenvalue.
    # And eigenvalues[1] is the largest eigenvalue.
    lambda_small = eigenvalues[0]
    lambda_large = eigenvalues[1]

    print("Step-by-step derivation:")
    print(f"1. The eigenvalues of the operator are calculated as:")
    print(f"   - Largest eigenvalue (lambda_large): {lambda_large:.4f}")
    print(f"   - Second-largest eigenvalue (lambda_small): {lambda_small:.4f}\n")

    print("2. The problem defines the post-measurement state as a superposition determined by the cubes of the eigenvalues.")
    print("3. We need to find the fidelity of this state with respect to the eigenstate of the second-largest eigenvalue (lambda_small).")
    print("4. The formula for this fidelity, F, simplifies to F = (lambda_small^6) / (lambda_large^6 + lambda_small^6).\n")
    
    # Step 3 & 4: Calculate the components of the simplified fidelity formula
    numerator = lambda_small**6
    denominator = lambda_large**6 + lambda_small**6
    
    fidelity = numerator / denominator

    # Step 5: Print the final equation with the computed numbers
    print("Final Equation Calculation:")
    print(f"The equation for fidelity is F = Numerator / Denominator")
    print(f"Numerator (lambda_small^6) = {numerator}")
    print(f"Denominator (lambda_large^6 + lambda_small^6) = {denominator}")
    
    print("\nResult:")
    print(f"Fidelity = {numerator} / {denominator}")
    print(f"The calculated fidelity is: {fidelity}")
    
    return fidelity

final_answer = solve_quantum_fidelity()
# The final answer is wrapped in <<<>>>
print(f"\n<<<{final_answer}>>>")
