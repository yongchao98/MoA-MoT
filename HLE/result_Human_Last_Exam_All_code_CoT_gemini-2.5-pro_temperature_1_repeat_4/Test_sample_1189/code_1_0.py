import math

# Please set the values for n and m here
n = 4
m = 4

class SetCounter:
    """
    A class to solve the specified combinatorial problem using a recurrence relation.
    """
    def __init__(self, n_val):
        """
        Initializes the solver for a given n.
        
        Args:
            n_val (int): The size of the base set S.
        """
        if not isinstance(n_val, int) or n_val <= 0:
            raise ValueError("n must be a positive integer.")
            
        self.n = n_val
        self.N = 2**n_val - 1
        # Memoization table for the recurrence function f(i)
        self.memo_f = {0: 1, 1: 0}

    def _f(self, i):
        """
        Calculates f(i), the number of valid sets of size i, using memoization.
        
        Args:
            i (int): The size of the set T.
            
        Returns:
            int: The value of f(i).
        """
        if i in self.memo_f:
            return self.memo_f[i]
        if i < 0:
            return 0

        # Recurrence relation:
        # i * f(i) = C(N, i-1) - f(i-1) - (N - i + 2) * f(i-2)
        
        # We use math.comb for efficient calculation of binomial coefficients.
        # It handles large numbers automatically.
        comb_val = math.comb(self.N, i - 1)
        
        f_i_minus_1 = self._f(i - 1)
        f_i_minus_2 = self._f(i - 2)
        
        numerator = comb_val - f_i_minus_1 - (self.N - i + 2) * f_i_minus_2
        
        # The result of counting sets must be an integer.
        result = numerator // i
        self.memo_f[i] = result
        return result

    def solve(self):
        """
        Solves the problem for the given n and m, and prints the detailed calculation.
        """
        if not isinstance(m, int) or m <= 0:
            raise ValueError("m must be a positive integer.")
        
        if m > self.N:
            print(f"For n={self.n}, m={m}:")
            print(f"The number of non-empty subsets is N = 2^{self.n}-1 = {self.N}.")
            print(f"It's impossible to choose m > N distinct subsets, so the answer is 0.")
            print(0)
            return
            
        # Calculate the final result by ensuring all necessary values are memoized
        final_result = self._f(m)

        print(f"For n={n}, m={m}:")
        
        if m <= 1:
            print(f"The number of sets is {final_result}.")
            return

        # For m >= 2, we can show the last step of the recurrence
        comb_val = math.comb(self.N, m - 1)
        fm1 = self.memo_f[m - 1]
        fm2 = self.memo_f[m - 2]
        coeff = self.N - m + 2
        
        print(f"The number of sets is f(m), calculated using the recurrence relation:")
        print(f"f(i) = (C(2^n-1, i-1) - f(i-1) - (2^n-i+1) * f(i-2)) / i")
        print("\nFor the final calculation step where i = {m}:")
        print(f"f({m}) = (C({self.N}, {m-1}) - f({m-1}) - ({self.N} - {m} + 2) * f({m-2})) / {m}")
        print(f"f({m}) = ({comb_val} - {fm1} - {coeff} * {fm2}) / {m}")
        print(f"f({m}) = ({comb_val - fm1} - {coeff * fm2}) / {m}")
        print(f"f({m}) = {(comb_val - fm1) - (coeff * fm2)} / {m}")
        print(f"f({m}) = {final_result}")


if __name__ == "__main__":
    try:
        solver = SetCounter(n)
        solver.solve()
    except ValueError as e:
        print(e)