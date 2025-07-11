def find_infimum_c():
    """
    This script finds the infimum of c for which the given Markov chain is transient.
    The derivation is explained step-by-step.
    """
    print("### Step-by-Step Derivation ###")
    
    # Step 1: Explain the methodology
    print("\nStep 1: Define the criterion for transience.")
    print("A Markov chain on the positive integers is transient if it has a sufficient drift towards infinity.")
    print("We use Lamperti's criterion for a random walk with state-dependent drift mu_k and variance sigma_k^2.")
    print("The chain is transient if the limit L = lim_{k->inf} (2 * k * mu_k) / sigma_k^2 is greater than 1.")
    
    # Step 2: Calculate the drift mu_k
    print("\nStep 2: Calculate the drift (mean displacement) mu_k for large k.")
    print("mu_k = E[X_{n+1} - k | X_n = k]")
    print("mu_k = (-2)*P_{k,k-2} + (-1)*P_{k,k-1} + (1)*P_{k,k+1} + (2)*P_{k,k+2}")
    print("mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print("mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2")
    print("mu_k = 2c/k")
    
    # Step 3: Calculate the variance sigma_k^2
    print("\nStep 3: Calculate the variance of displacement sigma_k^2 for large k.")
    print("First, we compute the second moment of displacement, E[Delta^2].")
    print("E[Delta^2] = (-2)^2*P_{k,k-2} + (-1)^2*P_{k,k-1} + (1)^2*P_{k,k+1} + (2)^2*P_{k,k+2}")
    print("E[Delta^2] = 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k) + 4*(1/4)")
    print("E[Delta^2] = 1 + 1/4 - c/k + 1/4 + c/k + 1")
    print("E[Delta^2] = 5/2")
    
    print("\nThe variance sigma_k^2 = E[Delta^2] - (mu_k)^2.")
    print("sigma_k^2 = 5/2 - (2c/k)^2")

    # Step 4: Apply the criterion
    print("\nStep 4: Evaluate the limit L for the transience criterion.")
    print("L = lim_{k->inf} (2 * k * mu_k) / sigma_k^2")
    print("As k -> inf, the (2c/k)^2 term in sigma_k^2 goes to 0, so sigma_k^2 approaches 5/2.")
    print("L = (2 * k * (2c/k)) / (5/2)")
    print("L = (4c) / (5/2)")
    print("L = 8c/5")
    
    # Step 5: Solve for c
    print("\nStep 5: Set up and solve the inequality for transience.")
    print("The chain is transient if L > 1.")
    # Outputting the numbers in the final equation
    print("The inequality is: (8 * c) / 5 > 1")
    print("Solving for c:")
    print("8 * c > 5")
    print("c > 5/8")

    # Step 6: State the final answer
    print("\nStep 6: Determine the infimum.")
    infimum_num = 5
    infimum_den = 8
    print(f"The set of values for c where the chain is transient is ( {infimum_num}/{infimum_den}, infinity ).")
    print(f"The infimum (greatest lower bound) of this set is {infimum_num}/{infimum_den}.")
    print(f"The value in decimal form is {infimum_num / infimum_den}.")

find_infimum_c()