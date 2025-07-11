import math

class TitanFraction:
    """
    Represents a fraction in the Titan 5-bit architecture.
    Numerators and denominators are constrained to the range [0, 31].
    Arithmetic operations automatically simplify results to fit constraints.
    """
    def __init__(self, num, den=1, verbose=True):
        if not (isinstance(num, int) and isinstance(den, int)):
            raise TypeError("Numerator and denominator must be integers.")
        if not (0 <= num <= 31 and 1 <= den <= 31):
            raise ValueError(f"Invalid 5-bit fraction: {num}/{den}. Must be in [0, 31].")
        
        self.num = num
        self.den = den
        self.verbose = verbose

    @staticmethod
    def _find_best_approx(val):
        """Finds the best 5-bit fraction approximation for a given decimal value."""
        if val == 0:
            return (0, 1)
            
        min_error = float('inf')
        best_frac = (0, 1)

        for b in range(1, 32):
            # Check the two integers closest to the ideal numerator
            a_ideal = b * val
            for a_candidate in [math.floor(a_ideal), math.ceil(a_ideal)]:
                if 0 <= a_candidate <= 31:
                    a = a_candidate
                    current_error = abs(val - a / b)
                    if current_error < min_error:
                        min_error = current_error
                        best_frac = (a, b)
        
        # Handle cases where the value is too large to be represented
        if val > 31 and best_frac == (0,1):
             return (31,1)

        return best_frac

    def _operate(self, other, op_name, op_func):
        """Generic method to handle binary operations."""
        n1, d1 = self.num, self.den
        n2, d2 = other.num, other.den
        
        new_num, new_den = op_func(n1, d1, n2, d2)
        
        op_symbol = {'add': '+', 'sub': '-', 'mul': '*', 'div': '/'}[op_name]
        
        if 0 <= new_num <= 31 and 1 <= new_den <= 31:
            if self.verbose:
                print(f"  {self} {op_symbol} {other} = {new_num}/{new_den}")
            return TitanFraction(new_num, new_den, verbose=self.verbose)
        else:
            if new_den == 0: raise ZeroDivisionError
            val = new_num / new_den
            approx_n, approx_d = self._find_best_approx(val)
            if self.verbose:
                print(f"  {self} {op_symbol} {other} = {new_num}/{new_den} -> too large, must approximate.")
                print(f"    Value is ~{val:.4f}. Best 5-bit approximation is {approx_n}/{approx_d}.")
            return TitanFraction(approx_n, approx_d, verbose=self.verbose)

    def __add__(self, other):
        return self._operate(other, 'add', lambda n1, d1, n2, d2: (n1 * d2 + n2 * d1, d1 * d2))

    def __sub__(self, other):
        return self._operate(other, 'sub', lambda n1, d1, n2, d2: (n1 * d2 - n2 * d1, d1 * d2))
        
    def __mul__(self, other):
        return self._operate(other, 'mul', lambda n1, d1, n2, d2: (n1 * n2, d1 * d2))

    def __truediv__(self, other):
        return self._operate(other, 'div', lambda n1, d1, n2, d2: (n1 * d2, d1 * n2))
        
    def __str__(self):
        return f"{self.num}/{self.den}"

    def __repr__(self):
        return f"TitanFraction({self.num}, {self.den})"
        
    @property
    def value(self):
        return self.num / self.den

def solve():
    print("--- Titan Computer Simulation: The Curious Monkey and Coconut Problem ---")
    
    # Step 1: Define constants as TitanFractions
    print("\n[Step 1] Defining constants with 5-bit fractional approximations:")
    pi = TitanFraction(22, 7)    # Approx for 3.14159...
    g = TitanFraction(29, 3)     # Approx for g = 9.8 m/s^2
    C45 = TitanFraction(12, 17)  # Approx for cos(45) = sqrt(2)/2 = 0.707...
    two = TitanFraction(2, 1)
    twenty = TitanFraction(20, 1)
    one = TitanFraction(1, 1)

    print(f"  π ≈ {pi} ({pi.value:.4f})")
    print(f"  g ≈ {g} ({g.value:.4f})")
    print(f"  cos(45°) ≈ {C45} ({C45.value:.4f})")
    
    # Step 2: Calculate mass (m)
    # m = (4/3) * π * r³ * ρ. With r=0.005m and ρ=900000 kg/m³, this simplifies to m = (3/20) * π
    print("\n[Step 2] Calculating the mass of the rock 'm = (3/20) * π':")
    three_twentieths = TitanFraction(3, 20)
    m = three_twentieths * pi
    print(f"  Final mass approximation: m = {m} kg (~{m.value:.4f} kg)")

    # Step 3: Calculate the required force F = (2 * g * m) / C45
    print("\n[Step 3] Calculating required force 'F'. We will try two values based on approximations.")
    # This calculation led to F values of 13/1 and 25/2 in the analysis. We test both.
    F_cand1 = TitanFraction(13, 1, verbose=False)
    F_cand2 = TitanFraction(25, 2, verbose=False)
    print(f"  Analysis shows two likely candidate forces due to approximation choices:")
    print(f"  Candidate Force 1: F₁ = {F_cand1} ({F_cand1.value:.4f} N)")
    print(f"  Candidate Force 2: F₂ = {F_cand2} ({F_cand2.value:.4f} N)")

    # Step 4 & 5: For each candidate F, calculate y and evaluate
    print("\n[Step 4 & 5] Simulating trajectory for each candidate Force.")
    hit = False
    for i, F in enumerate([F_cand1, F_cand2]):
        print(f"\n--- Testing F = {F} ---")
        print(f"Calculating landing height 'y = 20 * (1 - (g * m) / (F * cos(45°)))'")
        
        # Calculate term by term, showing Titan's operations
        print("  (a) Calculate g * m:")
        gm = g * m
        print("  (b) Calculate F * cos(45°):")
        F_C45 = F * C45
        print("  (c) Calculate T = (g*m) / (F*cos(45°)):")
        T = gm / F_C45
        print("  (d) Calculate 1 - T:")
        one_minus_T = one - T
        print("  (e) Calculate y = 20 * (1 - T):")
        y = twenty * one_minus_T
        
        y_val = y.value
        print(f"\n  Final calculated landing height for F={F}: y = {y} (~{y_val:.4f} m)")
        
        if 9.9 <= y_val <= 10.1:
            print(f"  SUCCESS! Landing height {y_val:.4f} m is within the target range [9.9, 10.1] m.")
            hit = True
            abs_error = abs(y_val - 10.0)
            print(f"<<<Y[{abs_error:.3f}]>>>")
            break
        else:
            print(f"  MISS! Landing height {y_val:.4f} m is outside the target range [9.9, 10.1] m.")

    if not hit:
        print("\n--- Conclusion ---")
        print("Neither of the feasible forces calculated by the Titan system resulted in a hit.")
        print("There are no other representable forces between 12.5 N and 13.0 N for the Titan computer to try.")
        print("Therefore, it is not possible to hit the coconut with this system.")
        print("<<<N0>>>")

solve()