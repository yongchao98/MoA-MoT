def calculate_chain_quality(beta, p):
    """
    Calculates the expected chain quality under a selfish mining attack.

    Args:
        beta (float): The portion of mining power held by the adversary (0 <= beta < 0.5).
        p (float): The probability that honest miners choose the adversary's block in a tie.
    """
    print(f"Calculating chain quality for beta = {beta} and p = {p}\n")

    if not (0 <= beta < 0.5):
        print("The model is valid for 0 <= beta < 0.5.")
        print("For beta >= 0.5, the selfish miner's chain grows infinitely longer, and chain quality approaches 0.")
        print("<<<0>>>")
        return
    
    if not (0 <= p <= 1):
        print("The probability 'p' must be between 0 and 1.")
        return

    # The formula for Chain Quality (Q) is:
    # Q = ((1-β) * (1-2β) * (1 + 2β - pβ)) / (1 - β - β²)

    # Calculate each part of the equation with the given numbers
    term1 = 1 - beta
    term2 = 1 - 2 * beta
    term3 = 1 + 2 * beta - p * beta
    
    numerator = term1 * term2 * term3

    den_term1 = beta
    den_term2 = beta**2
    denominator = 1 - den_term1 - den_term2
    
    # Print the equation with the calculated values
    print("Formula: Q = ((1-β) * (1-2β) * (1 + 2β - pβ)) / (1 - β - β²)\n")
    print("Plugging in the values:")
    print(f"Q = ( (1-{beta}) * (1-2*{beta}) * (1 + 2*{beta} - {p}*{beta}) ) / ( 1 - {beta} - {beta}² )")
    print(f"Q = ( ({term1}) * ({term2}) * ({term3}) ) / ( 1 - {den_term1} - {den_term2} )")
    print(f"Q = ( {term1} * {term2} * {term3} ) / ( {denominator} )")
    print(f"Q = {numerator} / {denominator}\n")
    
    if denominator == 0:
        print("Error: Denominator is zero. This occurs at beta = (sqrt(5)-1)/2, which is outside the valid range.")
        return
        
    final_quality = numerator / denominator

    print("Final Result:")
    print(f"The expected chain quality is: {final_quality}")
    print(f"<<<{final_quality}>>>")

if __name__ == '__main__':
    # --- User Parameters ---
    # β: Adversary's portion of mining power. Must be less than 0.5.
    beta = 0.3
    
    # p: Probability honest miners choose the adversary's block in a tie.
    p = 0.5
    # --- End of Parameters ---
    
    calculate_chain_quality(beta, p)