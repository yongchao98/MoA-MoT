import math

def calculate_c_k(n, k):
    """
    Calculates the Shapley value c_k for player p_k in the given game.
    
    The Shapley value is calculated as the expected marginal contribution of the player
    over all permutations. This can be expressed as:
    c_k = 4k*E[sigma_P^3] + 6k^2*E[sigma_P^2] + 4k^3*E[sigma_P] + k^4
    where sigma_P is the sum of indices of players preceding player k.
    
    This function computes the exact expectations by considering all possible
    predecessor sets P for player k.
    """
    if not (isinstance(n, int) and n > 1):
        raise ValueError("n must be an integer greater than 1.")
    if not (isinstance(k, int) and 1 <= k <= n):
        raise ValueError("k must be an integer between 1 and n.")

    # The set of other players' indices
    others = [i for i in range(1, n + 1) if i != k]
    
    # Calculate expectations by averaging over all possible predecessor sets.
    # The set of predecessors P for player k can have size s from 0 to n-1,
    # with each size having a probability of 1/n.
    # For a given size s, any subset of 'others' of size s is equally likely.
    
    e_sigma_p_1 = 0.0
    e_sigma_p_2 = 0.0
    e_sigma_p_3 = 0.0

    # Iterate over all possible sizes of the predecessor set
    for s in range(n):
        # The number of predecessor sets of size s
        num_subsets = math.comb(n - 1, s)
        
        # If there are no subsets of this size, the conditional expectation is 0
        if num_subsets == 0:
            e_s_1 = 0
            e_s_2 = 0
            e_s_3 = 0
        else:
            # We need the sum of sigma_P^r over all subsets P of size s
            sum_sigma_1 = 0
            sum_sigma_2 = 0
            sum_sigma_3 = 0
            
            # Iterate through all subsets of 'others' of size s
            # To do this, we iterate through bitmasks of length n-1 with s bits set
            # This is complex, so we use a simpler but correct expectation formula
            # E_s[sigma_P] = s/(n-1) * sum(others)
            # E_s[sigma_P^2] = s/(n-1) * sum(j^2 for j in others) + s(s-1)/((n-1)(n-2)) * (sum(others)^2 - sum(j^2 for j in others))
            # etc.

            # We can directly calculate the average contribution for a fixed predecessor size s
            
            # W_m(k) = sum of j^m for j != k
            W1 = sum(others)
            W2 = sum(j**2 for j in others)
            W3 = sum(j**3 for j in others)

            # Conditional expectations E_s[sigma_P^r] = E[sigma_P^r | |P|=s]
            e_s_1 = 0
            e_s_2 = 0
            e_s_3 = 0

            if n > 1 and s > 0:
                e_s_1 = (s / (n - 1)) * W1
            if n > 2 and s > 1:
                e_s_2 = (s / (n - 1)) * W2 + (s * (s - 1) / ((n - 1) * (n - 2))) * (W1**2 - W2)
            elif n==2 and s==1: # Handles case n-2=0 in denominator
                 e_s_2 = (s / (n-1)) * W2
            elif s==1:
                 e_s_2 = (s/(n-1)) * W2
            
            if n > 3 and s > 2:
                term1 = (s / (n-1)) * W3
                term2 = (3 * s * (s-1) / ((n-1)*(n-2))) * (W1 * W2 - W3)
                term3_coeff = (s * (s-1) * (s-2) / ((n-1)*(n-2)*(n-3)))
                term3_val = W1**3 - 3*W1*W2 + 2*W3
                term3 = term3_coeff * term3_val
                e_s_3 = term1 + term2 + term3
            elif s==1:
                 e_s_3 = (s / (n-1)) * W3
            elif s==2:
                 e_s_3 = (s / (n-1)) * W3 + (3 * s * (s-1) / ((n-1)*(n-2))) * (W1*W2-W3)

        # The probability of a predecessor set of size s is 1/n
        prob_s = 1 / n
        
        e_sigma_p_1 += prob_s * e_s_1
        e_sigma_p_2 += prob_s * e_s_2
        e_sigma_p_3 += prob_s * e_s_3

    c_k = 4*k*e_sigma_p_3 + 6*k*k*e_sigma_p_2 + 4*k*k*k*e_sigma_p_1 + k*k*k*k
    return c_k

# Example for n=3, to be used in the final formula output
n_val = 3
# Player 1's share
k1 = 1
c1 = calculate_c_k(n_val, k1)
# Player 2's share
k2 = 2
c2 = calculate_c_k(n_val, k2)
# Player 3's share
k3 = 3
c3 = calculate_c_k(n_val, k3)

print(f"For n = {n_val}, the total earnings are {int(sum([c1, c2, c3]))} dollars.")
print("The fair division is as follows:")
print(f"p_1 gets c_1 = {int(c1)} dollars.")
print(f"p_2 gets c_2 = {int(c2)} dollars.")
print(f"p_3 gets c_3 = {int(c3)} dollars.")
print("\nThe general formula for c_k given a fixed n is implemented in the python code.")
print("Let S be a subset of players {p_1, ..., p_n} excluding p_k.")
print("Let sigma_S be the sum of indices of players in S (e.g., if S={p_1,p_3}, sigma_S=1+3=4).")
print("The formula for c_k is the sum over all such subsets S:")
print("c_k = SUM_S [ (|S|! * (n-|S|-1)!) / n! * ( (sigma_S + k)^4 - (sigma_S)^4 ) ]")
print("\nFor example, when n=3, p_2's share (c_2) is calculated as:")
# Show one term of the sum for n=3, k=2
# Let's take S={p_1}
S_card = 1
n = 3
k = 2
sigma_S = 1
weight = (math.factorial(S_card) * math.factorial(n - S_card - 1)) / math.factorial(n)
marginal_cont = (sigma_S + k)**4 - sigma_S**4
term_val = weight * marginal_cont
print(f"Term for S={{p_1}}: (1! * (3-1-1)!)/3! * ((1 + 2)^4 - 1^4) = {weight:.3f} * ({int((1+2)**4)} - {int(1**4)}) = {term_val:.2f}")
