def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The portion of mining power held by the adversary (0 < beta < 1).
        p (float): The probability that honest miners choose the adversary's block in a tie (0 <= p <= 1).
    """

    print(f"Inputs:")
    print(f"Adversary's mining power (β) = {beta}")
    print(f"Honest miners' choice probability (p) = {p}\n")

    if beta >= 0.5:
        print("Since the adversary's mining power (β) is 0.5 or greater, their lead over the honest chain will grow indefinitely.")
        print("In the long run, they will author all blocks on the longest chain.")
        chain_quality = 0
        print("\n--- Final Equation ---")
        print(f"Expected Chain Quality = {chain_quality}")
        print("\n<<<0>>>")
        return

    alpha = 1.0 - beta
    print(f"Step 1: Calculate honest mining power (α)")
    print(f"α = 1 - β = 1 - {beta} = {alpha}\n")

    # In our model, the steady state probabilities follow a geometric distribution
    # for states i >= 1, with ratio q = beta / alpha.
    q = beta / alpha
    print(f"Step 2: Calculate the ratio q for the geometric distribution of states")
    print(f"q = β / α = {beta} / {alpha} = {q}\n")

    # These rates are derived from the steady-state analysis of the underlying Markov chain model.
    # λ_h is the long-term rate of honest blocks being added to the main chain.
    # λ_a is the long-term rate of adversary blocks being added to the main chain.
    
    # λ_h = (α - β) * (1 + β * (2 - p))
    # λ_a = (1/α) * (β**2 + β * (α - β) * (2*β + α*p))

    # We print the calculation step-by-step
    print("Step 3: Calculate the rate of honest blocks added to the chain (λ_h)")
    term_h_1 = alpha - beta
    term_h_2 = 1 + beta * (2 - p)
    lambda_h = term_h_1 * term_h_2
    print(f"λ_h = (α - β) * (1 + β * (2 - p))")
    print(f"λ_h = ({alpha} - {beta}) * (1 + {beta} * (2 - {p}))")
    print(f"λ_h = {term_h_1} * {term_h_2} = {lambda_h}\n")
    
    print("Step 4: Calculate the rate of adversary blocks added to the chain (λ_a)")
    term_a_1 = beta**2
    term_a_2 = beta * (alpha - beta)
    term_a_3 = 2*beta + alpha*p
    lambda_a = (1/alpha) * (term_a_1 + term_a_2 * term_a_3)
    print(f"λ_a = (1/α) * [β^2 + β * (α - β) * (2*β + α*p)]")
    print(f"λ_a = (1/{alpha}) * [{beta}^2 + {beta} * ({alpha} - {beta}) * (2*{beta} + {alpha}*{p})]")
    print(f"λ_a = (1/{alpha}) * [{term_a_1} + {term_a_2} * {term_a_3}] = {lambda_a}\n")

    print("Step 5: Calculate the final Expected Chain Quality")
    print("Chain Quality = λ_h / (λ_h + λ_a)")
    total_lambda = lambda_h + lambda_a
    chain_quality = lambda_h / total_lambda
    print(f"Chain Quality = {lambda_h} / ({lambda_h} + {lambda_a})")
    print(f"Chain Quality = {lambda_h} / {total_lambda} = {chain_quality}\n")
    
    print("--- Final Equation ---")
    print(f"Expected Chain Quality = ({alpha} - {beta}) * (1 + {beta} * (2 - {p})) / ( (({alpha} - {beta}) * (1 + {beta} * (2 - {p}))) + ((1/{alpha}) * ({beta**2} + ({beta} * ({alpha} - {beta})) * (2*{beta} + {alpha}*{p}))) ) = {chain_quality}")
    
    print(f"\n<<<{chain_quality}>>>")


if __name__ == '__main__':
    # --- User-defined variables ---
    # Adversary's mining power (portion of total). Must be < 1.
    BETA = 0.3
    # Probability honest miners choose the adversary's block in a tie.
    P = 0.5
    # --- End of user-defined variables ---

    if not (0 < BETA < 1 and 0 <= P <= 1):
        print("Error: Please set BETA to be between 0 and 1, and P to be between 0 and 1.")
    else:
        calculate_chain_quality(BETA, P)