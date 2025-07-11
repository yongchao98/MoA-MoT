import math

def get_transition_probs(k):
    """Calculates transition probabilities for state k >= 1."""
    if k == 0:
        return 1.0, 0, 0 # p(0,1)=1
    p_k = 0.5 * ((k + 1) / (k + 2))**3
    q_k = 0.5
    r_k = 1 - p_k - q_k
    return p_k, q_k, r_k

def f(k):
    """Calculates the value of the function f at state k."""
    if k == 0:
        return 0
    # Sum of (j+1)^3/64 for j=1 to k
    # which is 1/64 * ( ( (k+1)(k+2)/2 )**2 - 1 )
    term = (k + 1) * (k + 2) / 2
    return (term**2 - 1) / 64.0

def verify_counterexample(k):
    """
    Verifies that for a given state k, the drift condition holds.
    The drift is E[f(X_1)|X_0=k] - f(k).
    """
    if k == 0:
        print("k=0 is in set A, condition does not need to be checked.")
        return

    print(f"--- Verification for state k = {k} ---")
    
    p_k, q_k, r_k = get_transition_probs(k)
    f_k_minus_1 = f(k - 1)
    f_k = f(k)
    f_k_plus_1 = f(k + 1)

    print(f"For k={k}:")
    print(f"  p(k,k+1) = p_{k} = {p_k:.6f}")
    print(f"  p(k,k-1) = q_{k} = {q_k:.6f}")
    print(f"  p(k,k)   = r_{k} = {r_k:.6f}\n")
    
    print("Function values:")
    print(f"  f(k-1) = f({k-1}) = {f_k_minus_1:.6f}")
    print(f"  f(k)   = f({k}) = {f_k:.6f}")
    print(f"  f(k+1) = f({k+1}) = {f_k_plus_1:.6f}\n")

    # Calculate the drift
    drift = p_k * f_k_plus_1 + q_k * f_k_minus_1 + r_k * f_k - f_k

    print("Checking the condition Sum[p(k,y)f(y)] - f(k) >= 0:")
    print(f"The equation is: (p_{k} * f(k+1)) + (q_{k} * f(k-1)) + (r_{k} * f(k)) - f(k)")
    print(f"Plugging in the numbers: ({p_k:.6f} * {f_k_plus_1:.6f}) + ({q_k:.6f} * {f_k_minus_1:.6f}) + ({r_k:.6f} * {f_k:.6f}) - {f_k:.6f}")
    
    term1 = p_k * f_k_plus_1
    term2 = q_k * f_k_minus_1
    term3 = r_k * f_k
    print(f"Value of each term: ({term1:.6f}) + ({term2:.6f}) + ({term3:.6f}) - {f_k:.6f}")
    print(f"The result is: {drift:.6f}")

    if drift >= -1e-9: # Use a small tolerance for floating point errors
        print("\nThe condition is satisfied (result is >= 0).")
    else:
        print("\nThe condition is NOT satisfied.")

# Run verification for a specific state k, where k is not in A={0}.
verify_counterexample(k=5)
