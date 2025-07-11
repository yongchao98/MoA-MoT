import math

def main():
    """
    Calculates the number of different score sequences possible in a tournament
    between two teams.
    """
    
    print("The problem is to find the number of integer sequences (a_1..a_8, b_1..b_6) satisfying:")
    print("1. 0 <= a_i <= 6  (Team A has 8 players, each plays 6 games)")
    print("2. 0 <= b_j <= 8  (Team B has 6 players, each plays 8 games)")
    print("3. Sum(a_i) + Sum(b_j) = 48 (Total games played)")
    
    print("\nThis problem can be solved using generating functions. The total number of sequences")
    print("is the coefficient of x^48 in the expansion of P(x) = (1+...+x^6)^8 * (1+...+x^8)^6.")
    
    print("\nThe coefficient is calculated by the formula:")
    print("Total = Sum_{i=0 to 8, j=0 to 6} (-1)^(i+j) * C(8, i) * C(6, j) * C(61 - 7*i - 9*j, 13)")
    print("\nCalculating the individual non-zero terms of this sum:")
    
    def nCr(n, r):
        """
        Calculates the binomial coefficient C(n, r). Returns 0 if r < 0 or r > n.
        """
        if r < 0 or r > n:
            return 0
        return math.comb(n, r)

    total_sequences = 0
    
    # Iterate over the indices i and j from the formula
    for i in range(9):  # Corresponds to the (1-x^7)^8 term, so i from 0 to 8
        for j in range(7):  # Corresponds to the (1-x^9)^6 term, so j from 0 to 6
            
            # The 'n' for the final C(n, k) term, where k is 13
            n = 61 - 7*i - 9*j
            k = 13
            
            # A term is non-zero only if n >= k
            if n >= k:
                sign = (-1)**(i + j)
                c1 = nCr(8, i)
                c2 = nCr(6, j)
                c3 = nCr(n, k)
                
                term = sign * c1 * c2 * c3
                
                # We only care about non-zero terms
                if term != 0:
                    sign_char = "+" if sign > 0 else "-"
                    print(f"Term for (i={i}, j={j}): {sign_char} C(8,{i})*C(6,{j})*C({n},{k}) = {term}")
                    total_sequences += term
                    
    print("\n-------------------------------------------------------------")
    print(f"The final sum gives the total number of different sequences.")
    print(f"Total sequences = {total_sequences}")
    print("-------------------------------------------------------------")

if __name__ == "__main__":
    main()
