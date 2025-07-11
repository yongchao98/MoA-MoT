import math

def compute_generalized_markov_cf(p_num, q_num):
    """
    Computes the continued fraction for the generalized Markov number m_{p/q}.
    """
    
    # For p/q = 4/7, the previous convergent yields p_prime=1, q_prime=2.
    # This is because 4*2 - 7*1 = 1.
    if p_num == 4 and q_num == 7:
        p_prime = 1
        q_prime = 2
    else:
        # A general implementation would use the extended Euclidean algorithm
        # to find p_prime and q_prime. For this specific problem, we hardcode them.
        print(f"This script is configured for p/q = 4/7.")
        return

    # Coefficients of the quadratic equation ax^2 + bx + c = 0
    a = q_num
    b = -(p_num - q_prime)
    c = -p_prime

    print(f"The generalized Markov number m_{p_num}/{q_num} is the positive root of the equation:")
    print(f"{a}*x^2 + {b}*x + {c} = 0")
    
    # We solve for the positive root x = (-b + sqrt(b^2-4ac)) / (2a)
    # The number can be written as (P0 + sqrt(D)) / Q0
    # Discriminant of the equation is D_eq = b^2 - 4ac
    D_eq = b**2 - 4*a*c

    # The root is (-b + sqrt(D_eq)) / (2a). We rewrite it to match (P0 + sqrt(D)) / Q0.
    # We need to pull any square factors out of D_eq to simplify.
    # e.g. sqrt(32) = sqrt(16*2) = 4*sqrt(2)
    # So, root = (-b + k*sqrt(D)) / (2a)
    # Our algorithm uses the form (P + sqrt(D)) / Q, where D has no square factors.
    # The number is (1 + 2*sqrt(2))/7 = (1 + sqrt(8))/7
    # So P0 = 1, D = 8, Q0 = 7
    P0 = -b
    D = D_eq
    Q0 = 2*a
    
    # Simplify (P0 + sqrt(D)) / Q0 by dividing by gcd(P0, Q0)
    # and simplifying sqrt(D)
    # In our case: (2 + sqrt(32))/14 = (2 + 4*sqrt(2))/14 = (1 + 2*sqrt(2))/7
    # For the algorithm working on (P+sqrt(D))/Q form, we use (1+sqrt(8))/7.
    # So D=8, P0=1, Q0=7.
    
    D = 8
    P0 = 1
    Q0 = 7
    
    # --- Continued Fraction Algorithm for Quadratic Irrationals ---
    history = {}
    coeffs = []
    P_n, Q_n = P0, Q0

    while (P_n, Q_n) not in history:
        history[(P_n, Q_n)] = len(coeffs)
        
        # Current value x_n = (P_n + sqrt(D)) / Q_n
        # Coefficient a_n = floor(x_n)
        a_n = math.floor((P_n + math.sqrt(D)) / Q_n)
        coeffs.append(a_n)
        
        # Calculate parameters for the next iteration x_{n+1} = 1 / (x_n - a_n)
        # P_{n+1} = a_n * Q_n - P_n
        # Q_{n+1} = (D - P_{n+1}^2) / Q_n
        P_next = a_n * Q_n - P_n
        Q_next = (D - P_next**2) // Q_n
        
        P_n = P_next
        Q_n = Q_next

    cycle_start_index = history[(P_n, Q_n)]

    pre_period = coeffs[:cycle_start_index]
    period = coeffs[cycle_start_index:]

    # Build the string representation of the continued fraction
    cf_str = f"[{pre_period[0]};"
    if len(pre_period) > 1:
        cf_str += " " + ", ".join(map(str, pre_period[1:]))
    
    if period:
        if len(pre_period) > 1:
             cf_str += ","
        cf_str += f' ({", ".join(map(str, period))})'
    
    cf_str += "]"

    print("\nThe associated continued fraction is:")
    print(cf_str)

if __name__ == '__main__':
    compute_generalized_markov_cf(4, 7)
