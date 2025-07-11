import numpy as np
from scipy.optimize import fsolve

def find_valid_R_range():
    """
    This function calculates the valid range for the radius R of the set M
    for the Banach fixpoint theorem to apply.
    """
    
    # The two conditions for R are:
    # 1. exp(R) / 8 < 1  => R < ln(8) (Contraction)
    # 2. exp(R) / 8 <= R => exp(R) - 8R <= 0 (Invariance)

    # Condition 1: R must be less than ln(8)
    R_max = np.log(8)

    # Condition 2: Find the root of the equation exp(R) - 8R = 0
    # This equation has two roots. We are interested in the smaller positive root.
    # We define the function f(R) = exp(R) - 8R
    def f(R):
        return np.exp(R) - 8 * R

    # Use fsolve to find the root. We start the search near 0.1.
    R_min_root = fsolve(f, 0.1)[0]
    
    print("To prove the existence and uniqueness of the solution using the Banach fixpoint theorem,")
    print("we need to define a complete metric space M on which our operator T is a contraction.")
    print("\nLet's define M as a closed ball in the space of continuous functions satisfying the boundary conditions:")
    print("M = {u in C_0[0,1] | ||u||_inf <= R}")
    print("\nFor this to work, the radius R must satisfy two conditions:")
    print(f"1. Contraction: R < ln(8)  (approximately R < {R_max:.4f})")
    print(f"2. Invariance:   R >= R_0, where R_0 is the smallest positive root of exp(R) = 8R (approximately R >= {R_min_root:.4f})")

    print(f"\nThus, any R chosen from the interval [{R_min_root:.4f}, {R_max:.4f}) is a valid choice.")
    
    # We can choose any convenient number in this range. R=1 is a simple choice.
    chosen_R = 1.0
    is_valid = R_min_root <= chosen_R < R_max
    
    print(f"\nLet's pick a simple value for R, for example, R = {chosen_R}.")
    print(f"Is R={chosen_R} a valid choice? {is_valid}")
    
    print("\nTherefore, a suitable set M can be defined as:")
    
    # The final definition of M using the chosen R
    final_set_M = f"M = {{u ∈ C([0, 1]) | u(0) = u(1) = 0 and sup_{{x∈[0,1]}} |u(x)| ≤ {chosen_R}}}"
    print(final_set_M)


if __name__ == "__main__":
    find_valid_R_range()
