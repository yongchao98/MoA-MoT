def analyze_return_condition(p, p_prime):
    """
    Analyzes the 15s return condition based on a trade price (p)
    and the mid-price 15s later (p_prime).

    The problem states that the buyer's return is usually higher than the seller's return.
    This function demonstrates the mathematical implication of that statement.
    """
    
    # Definitions from the problem
    buyer_return = (p_prime / p) - 1
    seller_return = -((p_prime / p) - 1)
    
    print(f"Scenario analysis with trade price p = {p} and future mid-price p' = {p_prime}")
    print(f"Buyer's 15s Return: (p'/p - 1) = ({p_prime}/{p} - 1) = {buyer_return:.4f}")
    print(f"Seller's 15s Return: -(p'/p - 1) = -({p_prime}/{p} - 1) = {seller_return:.4f}")

    condition_met = buyer_return > seller_return
    print(f"Is Buyer's Return > Seller's Return? {condition_met}")
    
    # The inequality buyer_return > seller_return simplifies to p_prime > p
    # We print the final simplified equation (as an inequality)
    final_equation = "p' > p"
    print("\nThe observation that the buyer's return is higher than the seller's return is mathematically equivalent to the condition:")
    print(f"'{final_equation}'")
    print(f"Is this condition met in our scenario? {p_prime > p}\n")


# Case 1: The future price is higher
analyze_return_condition(p=100.0, p_prime=100.5)

# Case 2: The future price is lower
analyze_return_condition(p=100.0, p_prime=99.5)