def solve_chain_quality():
    """
    Calculates the expected chain quality under a specific selfish mining model.

    The function prompts the user for the adversary's mining power (beta)
    and the tie-breaking probability (p), then computes and prints the
    expected number of honest blocks, total blocks, and the final chain quality.
    """
    try:
        beta_str = input("Enter the adversary's mining power portion (β, e.g., 0.3): ")
        beta = float(beta_str)
        if not (0 <= beta < 1):
            print("Error: β must be a value between 0 (inclusive) and 1 (exclusive).")
            return

        p_str = input("Enter the tie-breaking probability (p, e.g., 0.5): ")
        p = float(p_str)
        if not (0 <= p <= 1):
            print("Error: p must be a value between 0 and 1 (inclusive).")
            return

    except ValueError:
        print("Invalid input. Please enter numerical values.")
        return

    # If beta is 0, the mining is entirely honest.
    if beta == 0:
        expected_honest_blocks = 1.0
        expected_total_blocks = 1.0
        chain_quality = 1.0
    else:
        alpha = 1.0 - beta

        # E[H] = α * (1 + 2αβ - αβp)
        expected_honest_blocks = alpha * (1 + 2 * alpha * beta - alpha * beta * p)
        
        # E[Total] = (1 - β² + β³) / α
        expected_total_blocks = (1 - beta**2 + beta**3) / alpha

        if expected_total_blocks == 0:
            chain_quality = 0.0
        else:
            chain_quality = expected_honest_blocks / expected_total_blocks

    print("\n--- Calculation Results ---")
    print(f"For β = {beta} and p = {p}:")
    print("\nThe chain quality is the ratio of expected honest blocks to expected total blocks added to the chain per 'round'.")
    print("\nChain Quality = (Expected Honest Blocks) / (Expected Total Blocks)")
    
    print("\nCalculated values for the equation:")
    print(f"Expected Honest Blocks = {expected_honest_blocks:.6f}")
    print(f"Expected Total Blocks  = {expected_total_blocks:.6f}")
    
    print(f"\nFinal Equation:")
    print(f"Expected Chain Quality = {expected_honest_blocks:.6f} / {expected_total_blocks:.6f}")
    
    print("\nResult:")
    print(f"The expected chain quality is: {chain_quality:.6f}")
    
    # Returning the final answer in the required format
    return f"<<<{chain_quality:.6f}>>>"

# Run the solver
# The final answer is wrapped in <<<>>> after the execution.
final_answer = solve_chain_quality()
if final_answer:
    print(final_answer)
