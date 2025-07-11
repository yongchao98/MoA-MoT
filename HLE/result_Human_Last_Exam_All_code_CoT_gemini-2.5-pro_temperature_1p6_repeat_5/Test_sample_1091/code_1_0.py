import math

def calculate_n_P_n(k):
    """
    Calculates n * P(n) for n = 6k.
    P(n) is the probability that the sum S is the zero vector.
    P(n) = (C(2k, k) / 2^(2k))^3
    n * P(n) = 6k * (C(2k, k) / 4^k)^3
    """
    n = 6 * k
    # Probability for one sum of 2k Rademacher variables to be 0
    # This uses math.comb for combinations C(n, k)
    try:
        # C(2k, k)
        comb = math.comb(2 * k, k)
        # P(X=0) = C(2k, k) / 2^(2k)
        p_zero = comb / (4**k)
        # P(n) = P(S=0) = P(X=0)^3
        prob_n = p_zero**3
        # n * P(n)
        result = n * prob_n
        return result, comb
    except (ValueError, OverflowError):
        # For very large k, use log gamma for better precision
        # log(P(X=0)) = log(C(2k,k)) - 2k*log(2)
        # log(C(2k,k)) = log( (2k)! ) - 2*log( k! )
        # log(m!) = gammaln(m+1)
        log_comb = math.lgamma(2 * k + 1) - 2 * math.lgamma(k + 1)
        log_p_zero = log_comb - (2 * k) * math.log(2)
        log_prob_n = 3 * log_p_zero
        log_result = math.log(n) + log_prob_n
        return math.exp(log_result), None

def main():
    """
    Main function to demonstrate the limit calculation.
    """
    print("This script computes n*P(n) for n=6k, where P(n) is the probability that ||S|| <= sqrt(2).")
    print("As shown in the analysis, this is equivalent to calculating n*P(S=0).")
    print("The formula is: n*P(n) = 6k * (C(2k, k) / 2^(2k))^3\n")
    
    print("The limit of n*P(n) as n (and k) go to infinity is 0.")
    print("Let's calculate the value for increasing k to observe the trend:\n")

    print(f"{'k':>5} {'n':>8} {'C(2k, k)':>20} {'n*P(n)':>25}")
    print("="*60)
    
    for k in [1, 2, 5, 10, 20, 50, 100, 200]:
        n = 6 * k
        result, comb_val = calculate_n_P_n(k)
        comb_str = str(comb_val) if comb_val is not None else "Overflow"
        print(f"{k:>5} {n:>8} {comb_str:>20} {result:>25.10f}")

if __name__ == "__main__":
    main()
