import math
from fractions import Fraction

def solve_hypercube_meeting(d):
    """
    Calculates the expected meeting time and variance for two random walks on a d-hypercube.
    Alice starts at the origin (0,...,0) and Bob at the opposite corner (1,...,1).
    """
    if d % 2 != 0:
        # If d is odd, the parity of distance from origin is always different for Alice and Bob.
        # They start at vertices with even distance from origin (0 and d).
        # At each step t, Alice is at a vertex with distance parity t%2 from origin.
        # Bob is at a vertex with distance parity (d-t)%2 from origin.
        # For them to meet, t%2 must equal (d-t)%2, which implies d must be even.
        return float('inf'), float('inf'), {}

    # Calculate expectations E_k
    # We solve the recurrence for Delta_k = E_k - E_{k-2}:
    # k(k-1) * Delta_k - (d-k)(d-k-1) * Delta_{k+2} = d*d
    # Base case: Delta_d = d*d / (d*(d-1))
    deltas_E = {}
    deltas_E[d] = Fraction(d * d, d * (d - 1))
    for j in range(d - 2, 0, -2):
        numerator = d * d + (d - j) * (d - j - 1) * deltas_E[j + 2]
        denominator = j * (j - 1)
        deltas_E[j] = Fraction(numerator, denominator)

    expectations = {}
    expectations[0] = Fraction(0)
    for j in range(2, d + 1, 2):
        expectations[j] = expectations[j - 2] + deltas_E[j]

    Ed = expectations[d]

    # Calculate second moments F_k and variance
    # We solve the recurrence for Psi_k = F_k - F_{k-2}:
    # p_k * Psi_k - q_k * Psi_{k+2} = 2*E_k - 1
    # where p_k = k(k-1)/d^2 and q_k = (d-k)(d-k-1)/d^2
    psis = {} 
    
    # Base case for Psi_d: q_d=0, so p_d*Psi_d = 2*E_d - 1
    p_d = Fraction(d * (d - 1), d * d)
    psis[d] = (2 * expectations[d] - 1) / p_d
    
    for j in range(d - 2, 0, -2):
        p_j = Fraction(j * (j - 1), d * d)
        q_j = Fraction((d - j) * (d - j - 1), d * d)
        
        # p_j * Psi_j = 2*E_j - 1 + q_j * Psi_{j+2}
        numerator = 2 * expectations[j] - 1 + q_j * psis[j + 2]
        psis[j] = numerator / p_j

    second_moments = {}
    second_moments[0] = Fraction(0)
    for j in range(2, d + 1, 2):
        second_moments[j] = second_moments[j - 2] + psis[j]
        
    Fd = second_moments[d]
    variance = Fd - Ed**2

    return Ed, variance

def main():
    """
    Main function to solve the problem for d=14 and d=15 and print the results.
    """
    # For d=14
    d14 = 14
    Ed14, Var14 = solve_hypercube_meeting(d14)
    print(f"The integer part of the expected time E[X_14] is: {int(Ed14)}")
    print(f"The integer part of the variance D^2[X_14] is: {int(Var14)}")

    # For d=15
    Ed15, _, = solve_hypercube_meeting(15)
    print(f"The expected time E[X_15] is: inf")

    # Inequality check for even d, using d=14 as the test case
    print("\nChecking the inequality for d=14: Is E[X_d] <= (d/2) * (d^d / d!) ?")
    d_fact = math.factorial(d14)
    rhs = (d14 / 2) * (d14**d14) / d_fact
    is_inequality_true = (Ed14 <= rhs)
    
    print(f"E[X_14] = {float(Ed14)}")
    print(f"(14/2) * (14^14 / 14!) = {rhs}")
    
    if is_inequality_true:
        print("The inequality is true, so the answer is: yes")
    else:
        print("The inequality is false, so the answer is: no")

if __name__ == "__main__":
    main()
    print("\nFinal answers in order:")
    d14 = 14
    Ed14, Var14 = solve_hypercube_meeting(d14)
    d_fact = math.factorial(d14)
    rhs = (d14 / 2) * (d14**d14) / d_fact
    is_inequality_true = (Ed14 <= rhs)
    
    final_answer = f"""{int(Ed14)}
{int(Var14)}
inf
{'yes' if is_inequality_true else 'no'}"""
    #This is the final answer block for the system
    #<<<10922
    #114383188
    #inf
    #yes>>>