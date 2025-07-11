import numpy as np
from scipy.optimize import fsolve

def define_set_M():
    """
    This function explains the derivation and calculates the constraints on the radius R
    for the set M used in the Banach fixed-point theorem proof for the given BVP.
    """

    print("To apply the Banach fixed-point theorem, we first reformulate the BVP as a fixed-point problem u = T(u).")
    print("This leads to the operator (Tu)(x) = ∫[0,1] G(x,s)exp(u(s))ds.\n")
    
    print("We choose the set M to be a closed ball in the space of continuous functions C[0,1]:")
    print("M = {u ∈ C[0, 1] | ||u||_∞ ≤ R} for some radius R > 0.\n")

    print("For T to be a contraction mapping on M, two conditions must be met:")
    print("1. Invariance (T(M) ⊆ M): Requires ||Tu||_∞ ≤ R. This leads to the condition: exp(R) / 8 ≤ R")
    print("2. Contraction: Requires the Lipschitz constant k < 1. This leads to the condition: exp(R) / 8 < 1  or  R < ln(8)\n")

    print("We now find the range of R satisfying both conditions.")
    
    # Define the function whose root we want to find for the invariance boundary
    # func(R) = 0  =>  R - exp(R)/8 = 0  =>  R = exp(R)/8
    func = lambda R_val: R_val - np.exp(R_val) / 8
    
    # We find the smaller positive root. A good initial guess is 0.1
    R_small_initial_guess = 0.1
    R_small = fsolve(func, R_small_initial_guess)[0]

    ln8 = np.log(8)
    
    print(f"The invariance condition is satisfied when R is between the two roots of R = exp(R)/8.")
    print(f"The smaller root is numerically found to be R ≈ {R_small:.5f}.")
    print(f"The contraction condition is R < ln(8) ≈ {ln8:.5f}.")
    print(f"Combining both, we need to choose a constant R such that {R_small:.5f} ≤ R < {ln8:.5f}.\n")

    print("--- The required set M is defined as follows ---")
    print("M = {u ∈ C[0, 1] : ||u||_∞ ≤ R}")
    print("where C[0, 1] is the Banach space of continuous functions on [0, 1] with the supremum norm,")
    print(f"and R is a constant chosen from the interval [{R_small:.5f}, {ln8:.5f}).")

define_set_M()