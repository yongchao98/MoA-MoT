import sys

# For floating point comparisons
TOL = 1e-9

class Q_sqrt2_sqrt3:
    """
    Represents an element of the field K = Q(sqrt(2), sqrt(3)).
    An element is represented as a + b*sqrt(2) + c*sqrt(3) + d*sqrt(6).
    """
    def __init__(self, a, b, c, d):
        # Allow integer and float coefficients
        self.coeffs = [float(a), float(b), float(c), float(d)]

    def __str__(self):
        return (f"({self.coeffs[0]:.2f}) + ({self.coeffs[1]:.2f})√2 + "
                f"({self.coeffs[2]:.2f})√3 + ({self.coeffs[3]:.2f})√6")

    def __eq__(self, other):
        return all(abs(x - y) < TOL for x, y in zip(self.coeffs, other.coeffs))

    def __mul__(self, other):
        a, b, c, d = self.coeffs
        ap, bp, cp, dp = other.coeffs
        
        # Coeff of 1:
        c0 = a*ap + 2*b*bp + 3*c*cp + 6*d*dp
        # Coeff of sqrt(2):
        c1 = a*bp + b*ap + 3*c*dp + 3*d*cp
        # Coeff of sqrt(3):
        c2 = a*cp + c*ap + 2*b*dp + 2*d*bp
        # Coeff of sqrt(6):
        c3 = a*dp + d*ap + b*cp + c*bp
        
        return Q_sqrt2_sqrt3(c0, c1, c2, c3)

# Define the automorphisms of K/Q
def sigma_2(x): # sends sqrt(2) -> -sqrt(2)
    return Q_sqrt2_sqrt3(x.coeffs[0], -x.coeffs[1], x.coeffs[2], -x.coeffs[3])

def sigma_3(x): # sends sqrt(3) -> -sqrt(3)
    return Q_sqrt2_sqrt3(x.coeffs[0], x.coeffs[1], -x.coeffs[2], -x.coeffs[3])

def sigma_4(x): # sends sqrt(2) -> -sqrt(2), sqrt(3) -> -sqrt(3)
    return sigma_2(sigma_3(x))

def check_orders():
    """
    Checks the orders of the extensions of automorphisms from K to L.
    An extension phi of sigma has order 2 if sigma(c)*c = 1,
    and order 4 if sigma(c)*c = -1, where c^2 = sigma(delta)/delta.
    """
    one = Q_sqrt2_sqrt3(1, 0, 0, 0)
    minus_one = Q_sqrt2_sqrt3(-1, 0, 0, 0)
    
    # Case 1: sigma_2
    # c_2^2 = (3 - 2*sqrt(2)) = (sqrt(2)-1)^2
    c2 = Q_sqrt2_sqrt3(-1, 1, 0, 0) # This is sqrt(2) - 1
    result = sigma_2(c2) * c2
    print("Checking extensions of sigma_2 (sends sqrt(2) -> -sqrt(2))")
    print(f"c2 = {c2}")
    print(f"sigma_2(c2) * c2 = {result}")
    if result == minus_one:
        print("Result is -1. All extensions of sigma_2 have order 4.\n")
    elif result == one:
        print("Result is 1. An extension of sigma_2 has order 2.\n")
        sys.exit("Error in logic: should be -1")

    # Case 2: sigma_3
    # c_3^2 = (2 - sqrt(3)) = ((sqrt(6)-sqrt(2))/2)^2
    c3 = Q_sqrt2_sqrt3(0, -0.5, 0, 0.5) # This is (sqrt(6)-sqrt(2))/2
    result = sigma_3(c3) * c3
    print("Checking extensions of sigma_3 (sends sqrt(3) -> -sqrt(3))")
    print(f"c3 = {c3}")
    print(f"sigma_3(c3) * c3 = {result}")
    if result == minus_one:
        print("Result is -1. All extensions of sigma_3 have order 4.\n")
    elif result == one:
        print("Result is 1. An extension of sigma_3 has order 2.\n")
        sys.exit("Error in logic: should be -1")

    # Case 3: sigma_4
    # c_4^2 = (3-2sqrt(2))(2-sqrt(3))
    c4 = c2 * c3
    result = sigma_4(c4) * c4
    print("Checking extensions of sigma_4 (sends sqrt(2)->-sqrt(2), sqrt(3)->-sqrt(3))")
    print(f"c4 = {c4}")
    print(f"sigma_4(c4) * c4 = {result}")
    if result == minus_one:
        print("Result is -1. All extensions of sigma_4 have order 4.\n")
    elif result == one:
        print("Result is 1. An extension of sigma_4 has order 2.\n")
        sys.exit("Error in logic: should be -1")

    print("Conclusion: All 6 automorphisms in G \\ H have order 4.")
    print("The group G has 1 element of order 1, 1 element of order 2, and 6 elements of order 4.")
    print("This is the signature of the Quaternion group Q_8.")

if __name__ == '__main__':
    check_orders()