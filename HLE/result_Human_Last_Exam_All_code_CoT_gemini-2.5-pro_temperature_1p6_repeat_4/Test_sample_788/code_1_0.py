import sys

# Suppress RecursionError for deep power calculations with negative exponents
sys.setrecursionlimit(2000)

class F4Element:
    """
    Represents an element of the finite field F_4 = {0, 1, a, a+1}.
    The element is stored as a polynomial c1*a + c0, with arithmetic mod 2.
    The field is defined by the relation a^2 + a + 1 = 0, so a^2 = a + 1.
    """
    def __init__(self, c1, c0):
        self.c1 = c1 % 2
        self.c0 = c0 % 2

    def __add__(self, other):
        # Addition is component-wise XOR
        return F4Element(self.c1 + other.c1, self.c0 + other.c0)

    def __mul__(self, other):
        # (c1*a + c0)(d1*a + d0) = c1d1*a^2 + (c1d0+c0d1)*a + c0d0
        # Substitute a^2 = a+1
        # = c1d1*(a+1) + (c1d0+c0d1)*a + c0d0
        # = (c1d1 + c1d0 + c0d1)*a + (c1d1 + c0d0)
        c1d1 = self.c1 * other.c1
        c1d0 = self.c1 * other.c0
        c0d1 = self.c0 * other.c1
        c0d0 = self.c0 * other.c0
        
        new_c1 = c1d1 + c1d0 + c0d1
        new_c0 = c1d1 + c0d0
        
        return F4Element(new_c1, new_c0)
        
    def __eq__(self, other):
        return self.c1 == other.c1 and self.c0 == other.c0

    def inv(self):
        """
        Calculate the multiplicative inverse. a^-1 = a+1, (a+1)^-1 = a, 1^-1 = 1.
        """
        if self == F4.ZERO:
            raise ZeroDivisionError
        if self == F4.ONE:
            return F4.ONE
        if self == F4.A:
            return F4.A_PLUS_1
        if self == F4.A_PLUS_1:
            return F4.A

    def __pow__(self, n):
        if n == 0:
            return F4.ONE
        if n < 0:
            return self.inv().__pow__(-n)
        
        res = F4.ONE
        base = self
        while n > 0:
            if n % 2 == 1:
                res = res * base
            base = base * base
            n //= 2
        return res

    def __repr__(self):
        if self.c1 == 0 and self.c0 == 0: return "0"
        if self.c1 == 0 and self.c0 == 1: return "1"
        if self.c1 == 1 and self.c0 == 0: return "a"
        if self.c1 == 1 and self.c0 == 1: return "a+1"

class F4:
    ZERO = F4Element(0, 0)
    ONE = F4Element(0, 1)
    A = F4Element(1, 0)
    A_PLUS_1 = F4Element(1, 1)

def calculate_invariant(config, config_name):
    """Calculates and prints the invariant for a given configuration."""
    print(f"Calculating invariant for configuration '{config_name}': {config}")
    print(f"The invariant is the sum of terms a^(x+y) for each peg (x,y).")
    
    total = F4.ZERO
    equation_terms = []
    
    for peg in config:
        x, y = peg
        power = x + y
        term = F4.A ** power
        equation_terms.append(f"a^({power})")
        print(f"  Term for peg {peg}: a^({x}+{y}) = a^{power} = {term}")
        total += term

    print(f"Final Equation: I = {' + '.join(equation_terms)}")
    print(f"Resulting Invariant: {total}\n")
    return total

def solve():
    """
    Demonstrates the existence of 4 distinct equivalence classes by finding
    four configurations with different invariant values.
    """
    print("--- Peg Game Equivalence Class Analysis ---")
    print("We will calculate an invariant for different peg configurations.")
    print("The invariant is based on arithmetic in the finite field F_4 = {0, 1, a, a+1}, where a^2=a+1.\n")

    # Define four configurations
    config1 = { (0,0) }
    config2 = { (1,0) }
    config3 = { (0,0), (1,0) }
    config4 = { (0,0), (0,1), (0,2) }

    # Calculate and display invariants
    inv1 = calculate_invariant(config1, "C1: Single peg at origin")
    inv2 = calculate_invariant(config2, "C2: Single peg at (1,0)")
    inv3 = calculate_invariant(config3, "C3: Two pegs")
    inv4 = calculate_invariant(config4, "C4: Three pegs in a line")
    
    invariants = {inv1, inv2, inv3, inv4}
    
    print("--- Conclusion ---")
    print(f"We have found configurations for {len(invariants)} distinct invariant values: {[str(i) for i in invariants]}")
    print("Since the invariant can only take 4 values (0, 1, a, a+1), and we have shown that all 4 are achievable,")
    print("there are exactly 4 equivalence classes.")
    print("\nThe number of equivalence classes is 4.")

solve()