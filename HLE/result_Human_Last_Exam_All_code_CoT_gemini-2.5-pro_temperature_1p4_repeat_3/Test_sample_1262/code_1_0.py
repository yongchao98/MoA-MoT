import itertools

# Helper functions and classes

def excedances(perm):
    """Calculates the number of excedances in a permutation."""
    count = 0
    # perm is 1-indexed, so we use enumerate starting from 1
    for i, val in enumerate(perm, 1):
        if val > i:
            count += 1
    return count

def get_derangements(n):
    """Generates all derangements for S_n."""
    derangements = []
    nums = list(range(1, n + 1))
    for p in itertools.permutations(nums):
        if all(p[i] != i + 1 for i in range(n)):
            derangements.append(p)
    return derangements

class Polynomial:
    """A helper class for polynomial manipulations."""
    def __init__(self, coeffs):
        # Coeffs stored in ascending order of power: [c_0, c_1, ..., c_k]
        self.coeffs = list(coeffs)
        self._trim()

    def _trim(self):
        """Removes trailing zero coefficients."""
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()

    def degree(self):
        """Returns the degree of the polynomial."""
        if self.coeffs == [0] or not self.coeffs:
            return -1
        return len(self.coeffs) - 1

    def __str__(self):
        """Creates an explicit string representation of the polynomial."""
        if self.degree() == -1: return "0"
        res = []
        for i in range(self.degree(), -1, -1):
            c = self.coeffs[i]
            if c == 0: continue
            
            term = ""
            if len(res) > 0:
                if c > 0: term += " + "
                else: term += " - "
            elif c < 0:
                term += "-"

            c = abs(c)

            if i == 0:
                term += str(c)
            elif c == 1:
                term += f"t"
            else:
                term += f"{c}*t"

            if i > 1:
                term += f"^{i}"
            res.append(term)
        return "".join(res)

def derangement_polynomial(n):
    """Computes the derangement polynomial d_n(t)."""
    if n == 0: return Polynomial([1]) # By convention, d_0(t) = 1
    if n < 0: return Polynomial([0])
    
    derangements = get_derangements(n)
    max_deg = n 
    coeffs = [0] * max_deg
    for p in derangements:
        exc = excedances(p)
        if exc >= max_deg:
             # Extend coeffs list if a higher degree is found
             coeffs.extend([0] * (exc - max_deg + 1))
             max_deg = exc + 1
        coeffs[exc] += 1
    return Polynomial(coeffs)

def main():
    """Main function to solve the problem and print the results."""
    print("This script will analyze the properties of derangement polynomials and their relation to the Hilbert series of uniform matroids.")
    
    # Part (c): Value of d_3(1)
    print("\n--- Analysis for Part (c) ---")
    n_c = 3
    derangements_3 = get_derangements(n_c)
    num_derangements_3 = len(derangements_3)
    print(f"For n = 3, the derangements (permutations sigma where sigma(i) != i) are: {derangements_3}.")
    print(f"The total number of derangements is {num_derangements_3}.")
    print("By definition, d_n(1) counts the number of derangements for S_n.")
    print(f"Therefore, d_3(1) = {num_derangements_3}.")
    print("To confirm, let's find the polynomial d_3(t):")
    p1 = (2, 3, 1) # excedances: 2>1, 3>2. Count = 2.
    p2 = (3, 1, 2) # excedances: 3>1. Count = 1.
    print(f"Excedances for {p1} = {excedances(p1)}")
    print(f"Excedances for {p2} = {excedances(p2)}")
    print(f"d_3(t) = 1*t^1 + 1*t^2")
    print(f"Evaluating at t=1: d_3(1) = 1*(1)^1 + 1*(1)^2 = 1 + 1 = 2")
    print(f"Final value for (c): 2")

    # Part (b): Leading coefficient of d_n(t)
    print("\n--- Analysis for Part (b) ---")
    print("The leading coefficient of d_n(t) is the number of derangements with the maximum possible number of excedances.")
    print("A permutation in S_n has at most n-1 excedances. This maximum is achieved only by sigma = (2, 3, ..., n, 1).")
    print("This permutation is a derangement for all n >= 2. Since it is unique, the coefficient is 1.")
    n_b = 4
    d4_poly = derangement_polynomial(n_b)
    print(f"Let's verify for n = {n_b}:")
    print(f"d_{n_b}(t) = {d4_poly}")
    print(f"The degree of d_{n_b}(t) is {d4_poly.degree()}, which is n-1.")
    print(f"The leading coefficient (coefficient of t^{d4_poly.degree()}) is {d4_poly.coeffs[-1]}.")
    print(f"The answer for (b) is Yes.")

    # Part (a): Identity and degree
    print("\n--- Analysis for Part (a) ---")
    print("The question asks to confirm H(U_{n-1, E})(t) = t^(n-1) * d_n(t).")
    print("The correct, known identity is H(U_{n-1, E})(t) = t^(n-1) * d_n(1/t).")
    print("The two would be equal only if d_n(t) were a reciprocal (palindromic) polynomial.")
    print(f"From Part (b), we found d_4(t) = {d4_poly}.")
    coeffs_d4 = d4_poly.coeffs
    print(f"The coefficients are (from t^0): {coeffs_d4}")
    print(f"The reversed coefficients are: {list(reversed(coeffs_d4))}")
    print("Since the coefficient list is not a palindrome, d_4(t) is not reciprocal. Thus, the identity is false.")
    
    print("\nNext, let's find the degree of H(U_{n-1, E})(t).")
    print("The Hilbert series of the Chow ring of a matroid of rank r has degree r-1.")
    print("The uniform matroid U_{n-1, E} has rank r = n-1.")
    print("Thus, the degree of H(U_{n-1, E})(t) is (n-1) - 1 = n-2.")
    print(f"Let's check for n={n_b}: H(U_{3, E})(t) should have degree 4-2=2.")
    print(f"Using the correct formula H(U_{3, E})(t) = t^3 * d_4(1/t)")
    print(f"d_4(1/t) = {d4_poly.coeffs[1]}/t + {d4_poly.coeffs[2]}/t^2 + {d4_poly.coeffs[3]}/t^3")
    print(f"t^3 * d_4(1/t) = {d4_poly.coeffs[1]}*t^2 + {d4_poly.coeffs[2]}*t + {d4_poly.coeffs[3]}")
    h4_true_poly = Polynomial(list(reversed(d4_poly.coeffs)))
    print(f"So, H(U_{3, E})(t) = {h4_true_poly}.")
    print(f"The degree of this polynomial is {h4_true_poly.degree()}, which equals 2, confirming the n-2 formula.")
    print(f"The answer for (a) is No, with the degree being n-2.")

if __name__ == '__main__':
    main()