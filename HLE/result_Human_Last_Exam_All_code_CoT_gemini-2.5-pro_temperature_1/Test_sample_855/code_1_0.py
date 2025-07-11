import math

class Fraction:
    """A simple fraction class to hold numerator and denominator."""
    def __init__(self, n, d=1):
        if not isinstance(n, int) or not isinstance(d, int):
            raise TypeError("Numerator and denominator must be integers.")
        if d == 0:
            raise ZeroDivisionError
        self.n = n
        self.d = d

    def __repr__(self):
        return f"({self.n}/{self.d})"

    @property
    def value(self):
        return self.n / self.d

class TitanNumber:
    """
    Represents a number in Titan's architecture using a fraction and a base-10 exponent.
    The fraction's numerator and denominator are constrained to be <= 31.
    """
    MAX_VAL = 31

    def __init__(self, fraction, exponent):
        self.fraction = fraction
        self.exponent = exponent
        self.simplify()

    def simplify(self):
        """Simplifies the fraction and adjusts the exponent if n or d > MAX_VAL."""
        # First, simple GCD simplification
        common = math.gcd(self.fraction.n, self.fraction.d)
        self.fraction.n //= common
        self.fraction.d //= common

        # Handle overflow by adjusting the exponent
        val = self.fraction.value
        if self.fraction.n > self.MAX_VAL or self.fraction.d > self.MAX_VAL:
            if val == 0: return
            
            # Adjust exponent so the value is between 1 and 10
            exp_change = math.floor(math.log10(abs(val)))
            self.exponent += exp_change
            new_val = val / (10**exp_change)
            
            # Find the best fractional approximation for the new value
            best_n, best_d = self.find_best_fraction(new_val)
            self.fraction.n = best_n
            self.fraction.d = best_d

    @staticmethod
    def find_best_fraction(value):
        """Finds a good n/d approximation for a float value."""
        min_error = float('inf')
        best_frac = (1, 1)
        for d_i in range(1, TitanNumber.MAX_VAL + 1):
            n_i = int(round(value * d_i))
            if 0 < n_i <= TitanNumber.MAX_VAL:
                error = abs(value - n_i / d_i)
                if error < min_error:
                    min_error = error
                    best_frac = (n_i, d_i)
        return best_frac

    @property
    def value(self):
        return self.fraction.value * (10**self.exponent)

    def __mul__(self, other):
        n_new = self.fraction.n * other.fraction.n
        d_new = self.fraction.d * other.fraction.d
        exp_new = self.exponent + other.exponent
        return TitanNumber(Fraction(n_new, d_new), exp_new)

    def __truediv__(self, other):
        n_new = self.fraction.n * other.fraction.d
        d_new = self.fraction.d * other.fraction.n
        exp_new = self.exponent - other.exponent
        return TitanNumber(Fraction(n_new, d_new), exp_new)

    def __repr__(self):
        return f"{self.fraction} * 10^{self.exponent} (val: {self.value:.2e})"

def main():
    # 1. Constants and Initial Values as TitanNumbers
    # Probe mass: 50 kg -> 5 * 10^1
    m_p = TitanNumber(Fraction(5, 1), 1)
    # Initial velocity: 300 m/s -> 3 * 10^2
    v_i = TitanNumber(Fraction(3, 1), 2)
    # Distance: 5000 m -> 5 * 10^3
    d = TitanNumber(Fraction(5, 1), 3)
    # Constant 0.5
    half = TitanNumber(Fraction(1, 2), 0)

    # 2. Calculate the deceleration force (F_decel = 0.5 * m * v_i^2 / d)
    v_i_sq = v_i * v_i  # (9/1)*10^4
    m_v_sq = m_p * v_i_sq # (45/1)*10^5 -> simplified to (9/2)*10^6
    half_m_v_sq = half * m_v_sq # (9/4)*10^6
    F_decel = half_m_v_sq / d # (9/20)*10^3

    # 3. Calculate the gravitational force (F_gravity = m_p * g)
    # First, calculate g = G * (4/3) * pi * R * rho
    # Approximations: G=6.67e-11~(20/3)*10^-11, pi~3.125=25/8
    # R_pandora = 2000km = 2e6m, rho_shell = 300kg/m^3 = 3e2
    TN_G = TitanNumber(Fraction(20, 3), -11)
    TN_4_3 = TitanNumber(Fraction(4, 3), 0)
    TN_pi = TitanNumber(Fraction(25, 8), 0)
    TN_R = TitanNumber(Fraction(2, 1), 6)
    TN_rho = TitanNumber(Fraction(3, 1), 2)

    # g calculation
    vc = TN_4_3 * TN_pi # (100/24)*10^0 -> (25/6)*10^0
    g = TN_G * vc * TN_R * TN_rho
    
    # F_gravity calculation
    F_gravity = m_p * g

    # 4. Calculate total rocket force (F_rocket = F_decel + F_gravity)
    # To add, convert to float, sum, and find best TN representation for the sum
    f_decel_val = F_decel.value
    f_grav_val = F_gravity.value
    f_rocket_val = f_decel_val + f_grav_val
    
    # For printing, round values to nearest integer
    f_decel_print = int(round(f_decel_val))
    f_grav_print = int(round(f_grav_val))
    f_rocket_print = f_decel_print + f_grav_print # Use sum of rounded parts for consistency

    # 5. Print the final equation
    # The final equation is F_rocket = F_decel + F_gravity
    print(f"{f_rocket_print} = {f_decel_print} + {f_grav_print}")

if __name__ == "__main__":
    main()