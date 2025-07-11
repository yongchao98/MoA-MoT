import numpy as np

def run_simulation(n):
    """
    Simulates the process for a given n to find a value for tau.
    """
    # Define problem parameters
    p = n**(-1/2)
    threshold = 1 - n**(-1/2)
    
    current_sum = 0
    # Loop from j=1 to n
    for j in range(1, n + 1):
        # Determine if X_j is 0 or U_j
        if np.random.rand() < p:
            # X_j = U_j
            u_j = np.random.uniform(0, n**(-1/2))
            current_sum += u_j
        else:
            # X_j = 0
            pass
            
        # Check stopping conditions
        if current_sum >= threshold:
            return j
        if j >= n:
            return n
    return n # Should be captured by the j >= n check above

def estimate_expectation(n, num_simulations=10000):
    """
    Estimates E[tau] by running many simulations.
    """
    taus = [run_simulation(n) for _ in range(num_simulations)]
    return np.mean(taus)

def analyze_theta():
    """
    Analyzes the problem to find theta based on the derivation.
    This code explains the rigorous mathematical argument rather than relying on simulation.
    """
    print("Step-by-step derivation for theta:")
    print("Let T = 1 - n^(-1/2) be the threshold.")
    print("The stopping time tau is min{j >= 1 : S_j >= T or j >= n}.")
    print("\nStep 1: Express E[tau] in terms of probabilities.")
    print("E[tau] = sum_{j=0}^{n-1} P(tau > j).")
    print("For j < n, P(tau > j) = P(S_k < T for all k <= j).")
    print("Since X_i >= 0, S_j is non-decreasing, so P(tau > j) = P(S_j < T).")
    print("E[tau] = sum_{j=0}^{n-1} (1 - P(S_j >= T)) = n - sum_{j=1}^{n-1} P(S_j >= T).")
    print("Thus, n - E[tau] = sum_{j=1}^{n-1} P(S_j >= T).")

    print("\nStep 2: Bound P(S_j >= T) using Chebyshev's inequality.")
    print("P(S_j >= T) <= Var(S_j) / (T - E[S_j])^2.")
    
    print("Calculating moments of X_i:")
    print("E[X_i] = (1/2) * n^(-1/2) * n^(-1/2) = 1/(2n).")
    print("E[X_i^2] = E[U_i^2] * n^(-1/2) = ( (n^(-1/2))^2 / 3 ) * n^(-1/2) = 1/(3n^(3/2)).")
    print("Var(X_i) = E[X_i^2] - (E[X_i])^2 = 1/(3n^(3/2)) - 1/(4n^2).")
    
    print("Mean and Variance of S_j = sum_{i=1 to j} X_i:")
    print("E[S_j] = j * E[X_i] = j/(2n).")
    print("Var(S_j) = j * Var(X_i) = j * (1/(3n^(3/2)) - 1/(4n^2)).")
    
    print("\nStep 3: Analyze the sum to find the order of magnitude of n - E[tau].")
    print("The bound is P(S_j >= T) <= [j * (1/(3n^(3/2)) - ...)] / [(1 - n^(-1/2) - j/(2n))^2].")
    print("For large n, this is approximately j/(3n^(3/2)) / (1-j/(2n))^2.")
    print("The sum sum_{j=1}^{n-1} P(S_j >= T) can be approximated by an integral:")
    print("Let j = x*n. The sum is on the order of n * integral(n*x/(3n^(3/2)) / (1-x/2)^2 dx) from 0 to 1.")
    print("This simplifies to O(n^(1/2)) * C, where C is a constant from the definite integral.")
    print("So, we have derived that n - E[tau] <= C * n^(1/2).")

    print("\nStep 4: Determine the value of theta.")
    print("The inequality E[tau] >= n - c*n^theta must hold.")
    print("We have shown E[tau] >= n - C*n^(1/2).")
    print("This implies the inequality holds for theta = 1/2.")
    theta = 1/2
    print(f"Theta is {theta}, which can be written as {int(theta*8)}/8.")
    
    numerator = int(theta * 8)
    denominator = 8
    
    print(f"\nThe value for theta is 1/2. We can write this as a multiple of 1/8.")
    print(f"The final equation to show this is:")
    print(f"{theta} = {numerator}/{denominator}")


analyze_theta()
final_answer = 1/2
<<<1/2>>>