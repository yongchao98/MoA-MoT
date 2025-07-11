import math

def calculate_expression(k, q):
    """
    Calculates the value of the expression for a given ideal power k and characteristic q.
    
    The ideal is m^k, where m is the maximal ideal of the ring of integers of a local field.
    N(I) = q^k
    phi(I) = q^k - q^(k-1)
    The expression is (phi(I)/N(I)) * log(log(N(I)))
    """
    if q <= 1:
        raise ValueError("q must be greater than 1.")
    if k <= 0:
        raise ValueError("k must be a positive integer.")
    if math.log(q) * k <= 1:
        # log(log(N(I))) is undefined if N(I) <= 1
        return float('nan')
        
    term = (1 - 1/q) * math.log(k * math.log(q))
    return term

def main():
    """
    Main function to demonstrate the growth of the expression.
    """
    # Let's use q=2 as an example for the characteristic.
    q = 2 
    
    print(f"Let q = {q}. We calculate the value of the expression for increasing k:")
    print("k\tValue")
    for k in [10, 100, 1000, 10000, 100000]:
        val = calculate_expression(k, q)
        print(f"{k}\t{val:.4f}")
    
    print("\nAs k increases, the expression grows indefinitely.")
    print("Therefore, the limit is positive infinity.")
    print("Final answer expressed as a string: '+inf'")

if __name__ == "__main__":
    main()
