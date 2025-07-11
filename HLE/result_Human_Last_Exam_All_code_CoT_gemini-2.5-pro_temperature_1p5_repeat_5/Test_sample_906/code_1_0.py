def solve_steady_state_pi0():
    """
    This function provides a step-by-step derivation for the steady-state 
    probability pi_0 of the given birth-death Markov process.
    """
    print("--- Derivation of the steady-state probability pi_0 ---")
    
    print("\nStep 1: State the general formula for steady-state probabilities in a birth-death process.")
    print("The steady-state probability pi_k can be expressed in terms of pi_0 as:")
    print("pi_k = pi_0 * Product_{i=0 to k-1} (lambda_i / mu_{i+1})")
    
    print("\nStep 2: Substitute the given rates into the formula.")
    print("The given rates are:")
    print("  - Birth rate lambda_i = lambda / (i + 1)")
    print("  - Death rate mu_{i+1} = mu")
    print("Substituting these into the formula gives:")
    print("pi_k = pi_0 * Product_{i=0 to k-1} [ (lambda / (i + 1)) / mu ]")

    print("\nStep 3: Simplify the expression using rho = lambda / mu.")
    print("pi_k = pi_0 * Product_{i=0 to k-1} [ (lambda / mu) * (1 / (i + 1)) ]")
    print("Let rho = lambda / mu. The expression becomes:")
    print("pi_k = pi_0 * Product_{i=0 to k-1} [ rho / (i + 1) ]")
    print("Expanding the product:")
    print("pi_k = pi_0 * [rho/(0+1)] * [rho/(1+1)] * ... * [rho/((k-1)+1)]")
    print("pi_k = pi_0 * (rho/1) * (rho/2) * ... * (rho/k)")
    print("This simplifies to:")
    print("pi_k = pi_0 * (rho^k / k!)")

    print("\nStep 4: Apply the normalization condition.")
    print("The sum of all probabilities must equal 1:")
    print("Sum_{k=0 to infinity} pi_k = 1")
    print("Substituting the expression for pi_k:")
    print("Sum_{k=0 to infinity} [pi_0 * (rho^k / k!)] = 1")
    print("Factoring out pi_0 gives:")
    print("pi_0 * (Sum_{k=0 to infinity} [rho^k / k!]) = 1")

    print("\nStep 5: Recognize the sum as a Taylor series.")
    print("The sum Sum_{k=0 to infinity} [rho^k / k!] is the Taylor series for e^rho.")
    
    print("\nStep 6: State the final equation and solve for pi_0.")
    print("Substituting the Taylor series sum into the normalization equation gives our final equation:")
    print("pi_0 * e^rho = 1")
    print("\nSolving for pi_0, we get the final result:")
    print("pi_0 = 1 / e^rho = e^(-rho)")

# Execute the derivation
solve_steady_state_pi0()