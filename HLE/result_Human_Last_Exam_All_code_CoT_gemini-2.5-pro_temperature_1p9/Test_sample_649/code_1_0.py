import math

def display_fanos_lower_bound():
    """
    This function outlines the derivation and presents the final lower bound on the minimax risk
    R*_n using Fano's inequality.
    """

    print("Deriving the lower bound for the minimax risk R*_n using Fano's method.")
    print("--------------------------------------------------------------------------\n")
    print("Step 1: Relate the estimation risk to the probability of testing error.")
    print("We choose a set of M+1 hypotheses {P_0, ..., P_M} such that for j!=k:")
    print("  (a) The parameters are separated: rho(theta_j, theta_k) >= delta")
    print("  (b) The distributions are close: KL(P_j || P_k) <= K")
    print("\nIf an estimator is accurate, rho(theta_hat, theta_j) < delta/2, it can correctly identify j.")
    print("This implies that the probability of making a large estimation error is linked to the probability of making a testing error.")
    print("Let P_e be the average probability of error in hypothesis testing.\n")
    
    print("Step 2: Link expected loss to testing error.")
    print("Using Markov's inequality, we can show that the Bayes risk is bounded by:")
    print("  Risk >= Phi(delta / 2) * P_e\n")

    print("Step 3: Apply Fano's inequality.")
    print("Fano's inequality provides a lower bound for P_e based on the mutual information I(S; J) between the data S and the hypothesis index J.")
    print("  P_e >= 1 - (I(S; J) + log(2)) / log(M+1)")
    print("\nThe mutual information for n i.i.d. samples is bounded by the KL divergence:")
    print("  I(S; J) <= n * max_{j,k} KL(P_j || P_k) = n * K\n")
    
    print("Step 4: Combine the pieces to get the final bound.")
    print("Substituting the bounds from steps 2 and 3, we arrive at the minimax lower bound:\n")

    # Final equation is R*_n >= Phi(delta/2) * (1 - (n*K + log(2)) / log(M+1))
    # We will print it component by component as requested.
    
    print("The final lower bound on R^*_n is:\n")
    print("R^*_n >= ", end="")
    print("Phi(delta / 2)", end="")
    print(" * ", end="")
    print("(1 - (n * K + log(2)) / log(M + 1))")
    
    print("\n\nWhere:")
    print("  R^*_n: The minimax risk.")
    print("  Phi: The non-decreasing loss shaping function.")
    print("  delta: The minimum separation between parameters in the chosen hypothesis set.")
    print("  n: The number of i.i.d. samples.")
    print("  K: The maximum Kullback-Leibler (KL) divergence between any pair of chosen distributions.")
    print("  M+1: The number of hypotheses carefully chosen from the family of distributions.")
    print("  log: The natural logarithm.")

# Execute the function to display the result
display_fanos_lower_bound()