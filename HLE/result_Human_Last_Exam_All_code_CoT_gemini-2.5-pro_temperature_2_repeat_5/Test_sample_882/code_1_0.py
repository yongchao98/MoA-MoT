import math

def solve():
    """
    Solves for p in the mixed strategy Nash Equilibrium and calculates the final requested value.
    """
    N = 8
    M = 3 * N

    # The equilibrium condition U_D(p) = U_S(p) simplifies to F(p) = 0.
    # F(p) = (1-p/N)^M - p*((N-1)p/N)^M - (1-p)
    # We solve this equation for p using a numerical method.
    def F(p):
        if p < 0 or p > 1:
            return float('nan')
        
        term1 = (1 - p / N)**M
        term2 = p * (((N - 1) * p / N)**M)
        term3 = 1 - p
        return term1 - term2 - term3

    # We need to solve the original form to avoid issues at p=0,1.
    # The numerically stable form is U_S(p) - U_D(p) = 0.
    def h(p):
        # A(p) corresponds to the term (1 - p/N)^M
        A = (1 - p / N)**M
        # B(p) corresponds to the term ((N-1)p/N)^M
        B = ((N - 1) * p / N)**M
        
        # Payoff component from Spread strategy
        # handle p=1 case with limit
        if p == 1:
            payoff_s_core = M * (1/N - 1) * ((N-1)/N)**(M-1) + ((N-1)/N)**M
        else:
            payoff_s_core = (A - B) / (1 - p)
        
        # Payoff component from Discrete strategy
        # handle p=0 case with limit
        if p == 0:
            payoff_d_core = M / N
        else:
            payoff_d_core = (1 - A) / p
            
        return payoff_s_core - payoff_d_core

    # Bisection method to find the root of h(p) = 0
    # Analytical checks confirm the root for N=8 is between 0.9 and 1.0
    low = 0.9
    high = 1.0
    
    for _ in range(100):  # 100 iterations provide sufficient precision
        mid = (low + high) / 2
        if mid == low or mid == high: # Stop if precision limit is reached
            break
        
        # Use a very small epsilon to avoid evaluating exactly at 1.0
        if mid >= 1.0:
            mid = 1.0 - 1e-15

        if h(mid) < 0:
            low = mid
        else:
            high = mid
            
    p_solution = (low + high) / 2

    # The equation being solved is (1 - p/8)^24 - p * (7*p/8)^24 = 1 - p
    print("The equilibrium probability p is the solution to the equation:")
    print("(1 - p/8)^24 - p * ((7*p)/8)^24 - (1-p) = 0")
    print(f"\nCalculated value of p: {p_solution:.10f}")

    # Calculate the final result as requested by the problem
    value = 10000 * (1 - p_solution)
    final_answer = math.floor(value)

    print(f"Value of 10000 * (1 - p): {value:.10f}")
    print(f"The floor of this value is: {final_answer}")
    return final_answer

result = solve()
print(f"<<<{result}>>>")