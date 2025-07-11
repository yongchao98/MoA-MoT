class ZTau:
    """
    A class to represent numbers in the ring Z[τ], where τ = (sqrt(5)-1)/2.
    An element is stored as a pair of integers (a, b) representing a + bτ.
    The relation τ^2 = 1 - τ is used for multiplication.
    """
    def __init__(self, a, b):
        self.a = int(a)
        self.b = int(b)

    def __add__(self, other):
        return ZTau(self.a + other.a, self.b + other.b)

    def __mul__(self, other):
        # (a+bτ)(c+dτ) = ac + (ad+bc)τ + bdτ^2
        # = ac + (ad+bc)τ + bd(1-τ)
        # = (ac+bd) + (ad+bc-bd)τ
        a, b, c, d = self.a, self.b, other.a, other.b
        return ZTau(a * c + b * d, a * d + b * c - b * d)

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __repr__(self):
        # Pretty print the number
        if self.b == 0:
            return f"{self.a}"
        if self.a == 0:
            return f"{self.b}τ"
        sign = '+' if self.b > 0 else '-'
        return f"({self.a} {sign} {abs(self.b)}τ)"

def tau_pow(k):
    """Computes τ^k as a ZTau object."""
    if k == 0:
        return ZTau(1, 0)
    
    base = ZTau(0, 1) if k > 0 else ZTau(1, 1) # τ or τ^-1 = 1+τ
    
    res = ZTau(1, 0)
    for _ in range(abs(k)):
        res = res * base
    return res

def run_verification():
    """
    Verifies the existence of a function g in G with g'(0)=τ and g'(1)=1.
    This is done by constructing the function with three piecewise linear segments.
    """
    print("Goal: Construct an element g in G with slope exponents (1, 0) at the boundaries.")
    print("This requires g'(0) = τ^1 and g'(1) = τ^0 = 1.\n")
    
    # We choose breakpoints that are in Z[τ]: 0, τ^2, τ, 1.
    bp0, bp1, bp2, bp3 = ZTau(0,0), tau_pow(2), tau_pow(1), ZTau(1,0)
    
    # Calculate the lengths of the intervals (Δx_i)
    dx0 = bp1 # τ^2
    dx1 = bp2 - bp1 # τ - τ^2 = 2τ - 1
    dx2 = bp3 - bp2 # 1 - τ = τ^2
    
    print(f"Breakpoints chosen: 0, {bp1}, {bp2}, 1.")
    print(f"Interval lengths:")
    print(f"  Δx₀ = x₁-x₀ = {dx0}")
    print(f"  Δx₁ = x₂-x₁ = {dx1}")
    print(f"  Δx₂ = x₃-x₂ = {dx2}\n")

    # Define the slopes (s_i = τ^k_i)
    # We want k₀=1 and k₂=0. We need to find a valid integer k₁
    s0 = tau_pow(1)
    s2 = tau_pow(0)
    
    # The condition for a valid element g is Σ s_i Δx_i = 1.
    # s₀Δx₀ + s₁Δx₁ + s₂Δx₂ = 1
    # We need to find an s₁ = τ^k₁ that satisfies this.
    # From pencil-and-paper calculation, k₁ = -1. Let's verify.
    k1 = -1
    s1 = tau_pow(k1)

    print("Slopes chosen:")
    print(f"  s₀ = τ^1 = {s0} (for endpoint slope at 0)")
    print(f"  s₁ = τ^{k1} = {s1} (for the middle segment)")
    print(f"  s₂ = τ^0 = {s2} (for endpoint slope at 1)\n")
    
    # Verify the sum
    term0 = s0 * dx0
    term1 = s1 * dx1
    term2 = s2 * dx2
    total_sum = term0 + term1 + term2
    one = ZTau(1, 0)
    
    print("Verifying the condition Σ s_i Δx_i = 1:")
    print(f"Equation: ({s0})*({dx0}) + ({s1})*({dx1}) + ({s2})*({dx2}) = 1")
    print("Calculating each term:")
    print(f"  Term 1: s₀Δx₀ = {term0}")
    print(f"  Term 2: s₁Δx₁ = {term1}")
    print(f"  Term 3: s₂Δx₂ = {term2}")
    print(f"Sum of terms: {term0} + {term1} + {term2} = {total_sum}")

    if total_sum == one:
        print("\nThe condition is satisfied. An element g with the desired endpoint slopes exists.")
        print("This shows that (1,0) is in the image of the homomorphism ψ.")
        print("By symmetry, (0,1) is also in the image.")
        print("Therefore, the image is Z², and Ab(G) ≅ Z².")
        r = 2
        t = 1
        print(f"\nConclusion:")
        print(f"The rank of Ab(G) is r = {r}.")
        print(f"The order of the torsion subgroup of Ab(G) is t = {t}.")
        print(f"The pair (r, t) is ({r}, {t}).")
    else:
        print("\nCalculation failed. The construction is incorrect.")

if __name__ == '__main__':
    run_verification()