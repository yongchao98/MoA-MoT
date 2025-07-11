# The problem asks for the value of alpha in the expression d_n = Theta(n^alpha).
# My step-by-step derivation shows that the asymptotic growth rate of d_n is Theta(n^6).
# By comparing the two expressions, we can determine the value of alpha.

# Let the asymptotic growth be d_n = C * n^alpha for some constant C.
# Our derivation concluded that d_n is proportional to n^6.
# Therefore, alpha must be 6.

alpha = 6
print(f"The asymptotic growth rate of d_n is Theta(n^alpha), where alpha = {alpha}.")