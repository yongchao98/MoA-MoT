import math

def get_divisors(n):
    """Helper function to get all divisors of a number."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def sigma(n, k):
    """Calculates the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    return sum(d**k for d in get_divisors(n))

class QSeries:
    """A class to handle q-series operations."""
    def __init__(self, coeffs):
        self.coeffs = list(coeffs)

    def __getitem__(self, key):
        return self.coeffs[key]

    def __len__(self):
        return len(self.coeffs)
        
    def __str__(self):
        res = []
        for i, c in enumerate(self.coeffs):
            if c != 0:
                if i == 0:
                    res.append(str(c))
                elif i == 1:
                    res.append(f"{c}q")
                else:
                    res.append(f"{c}q^{i}")
        return " + ".join(res)

    def __add__(self, other):
        n = max(len(self), len(other))
        new_coeffs = [0] * n
        for i in range(n):
            c1 = self.coeffs[i] if i < len(self) else 0
            c2 = other.coeffs[i] if i < len(other) else 0
            new_coeffs[i] = c1 + c2
        return QSeries(new_coeffs)

    def __sub__(self, other):
        return self + (other * -1)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return QSeries([c * other for c in self.coeffs])
        
        n = len(self) + len(other) - 1
        new_coeffs = [0] * n
        for i in range(len(self)):
            for j in range(len(other)):
                new_coeffs[i + j] += self.coeffs[i] * other.coeffs[j]
        return QSeries(new_coeffs)

def get_E_series(k, num_coeffs):
    """Generates the q-expansion for the normalized Eisenstein series E_k."""
    coeffs = [0] * num_coeffs
    coeffs[0] = 1
    
    if k == 4:
        prefactor = 240
        power = 3
    elif k == 8:
        prefactor = 480
        power = 7
    else:
        raise ValueError("Only k=4 and k=8 are supported.")

    for n in range(1, num_coeffs):
        coeffs[n] = prefactor * sigma(n, power)
    return QSeries(coeffs)

def transform_q2(series):
    """Transforms a q-series f(q) to f(q^2)."""
    n = len(series.coeffs)
    new_coeffs = [0] * (2*n - 1)
    for i in range(n):
      if 2*i < len(new_coeffs):
        new_coeffs[2*i] = series.coeffs[i]
    return QSeries(new_coeffs)
    
def solve():
    """
    Main function to solve the problem.
    """
    PRECISION = 10 # Number of coefficients to compute (q^0 to q^9)

    # 1. Generate E_4 and F
    E4 = get_E_series(4, PRECISION)
    F = transform_q2(E4)
    
    # 2. Generate the basis g1, g2, g3
    g1 = E4 * E4           # E_4^2
    g2 = E4 * F            # E_4 * F
    g3 = F * F             # F^2

    # 3. Construct basis for the cusp forms
    # These must have constant term 0
    f_a = g1 - g3
    f_b = g2 - g3

    # 4. Find the unique linear combination with a2=0
    a2_fa = f_a[2]
    a2_fb = f_b[2]
    
    # We create h = a2(f_a)*f_b - a2(f_b)*f_a, so a2(h) is 0.
    h = (f_b * a2_fa) - (f_a * a2_fb)
    
    # 5. Normalize h to get the final cusp form f
    a1_h = h[1]
    if a1_h == 0:
      raise ValueError("Constructed form has no q^1 term, cannot normalize.")
    
    f = h * (1 / a1_h)
    
    # 6. Find first three non-zero coefficients and their sum
    non_zero_coeffs = []
    print("The unique normalized cusp form is f(z) = ")
    equation_parts = []
    
    for i in range(1, len(f.coeffs)):
        coeff = round(f.coeffs[i])
        if coeff != 0:
            non_zero_coeffs.append(coeff)
            if i == 1:
                part = "q" if coeff == 1 else f"{coeff}*q"
            else:
                sign = "+" if coeff > 0 else "-"
                abs_coeff = abs(coeff)
                if len(equation_parts) > 0: # Add sign for subsequent terms
                    if abs_coeff == 1:
                        part = f" {sign} q^{i}"
                    else:
                        part = f" {sign} {abs_coeff}*q^{i}"
                else: # First term special handling
                     if coeff == 1:
                         part = f"q^{i}"
                     elif coeff == -1:
                         part = f"-q^{i}"
                     else:
                         part = f"{coeff}*q^{i}"
            
            equation_parts.append(part)
            if len(non_zero_coeffs) == 3:
                break
    
    print("".join(equation_parts).strip() + " + ...")
    
    # Print the equation representing the sum
    sum_val = sum(non_zero_coeffs)
    print("\nThe sum of the first three non-zero coefficients is:")
    sum_eq_parts = [str(c) for c in non_zero_coeffs]
    print(f"{' + '.join(sum_eq_parts)} = {sum_val}")
    
    return sum_val

# Execute the solution
result = solve()
# Final answer block. The calculated answer is wrapped as requested.
# The calculation shows the first three non-zero coeffs are 1, 140, 1024. Their sum is 1+140+1024=1165.
# Let me recheck calculation 1008000/7200=140. 7372800/7200=1024. Looks correct.

final_answer = int(result)
print(f"\n<<<${final_answer}$>>>")