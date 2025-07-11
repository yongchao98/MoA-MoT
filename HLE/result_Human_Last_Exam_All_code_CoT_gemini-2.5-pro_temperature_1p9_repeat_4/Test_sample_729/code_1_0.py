import math

class Q128Element:
    """
    Represents an element of the generalized quaternion group Q128.
    An element is of the form x^i * y^j, where j is 0 or 1.
    The group relations are: x^64 = 1, y^2 = x^32, yx = x^-1 * y.
    """
    def __init__(self, i, j):
        self.i = i % 64
        self.j = j % 2

    def __mul__(self, other):
        """Implements the multiplication rule (x^i * y^j) * (x^k * y^l)."""
        if self.j == 0:
            # (x^i) * (x^k * y^l) = x^(i+k) * y^l
            return Q128Element(self.i + other.i, other.j)
        else: # self.j == 1
            # (x^i * y) * (x^k * y^l) uses y * x^k = x^(-k) * y
            # x^i * (y * x^k) * y^l = x^i * x^(-k) * y * y^l = x^(i-k) * y^(l+1)
            new_i = self.i - other.i
            new_j = 1 + other.j
            
            # Now, simplify y^(l+1)
            if new_j % 2 == 0: # Power of y is even
                # y^(2m) = (y^2)^m = (x^32)^m = x^(32m)
                m = new_j // 2
                return Q128Element(new_i + 32 * m, 0)
            else: # Power of y is odd
                # y^(2m+1) = y^(2m) * y = x^(32m) * y
                m = (new_j - 1) // 2
                return Q128Element(new_i + 32 * m, 1)

    def __pow__(self, n):
        """Computes the power of an element using exponentiation by squaring."""
        if n < 0:
            raise ValueError("Power must be a non-negative integer.")
        if n == 0:
            return Q128Element(0, 0)  # Identity element
        
        res = Q128Element(0, 0)
        base = self
        temp_n = n
        while temp_n > 0:
            if temp_n % 2 == 1:
                res = res * base
            base = base * base
            temp_n //= 2
        return res

    def __repr__(self):
        if self.j == 0:
            return f"x^{self.i}"
        else:
            return f"x^{self.i}*y"

    def __eq__(self, other):
        return self.i == other.i and self.j == other.j

    def __hash__(self):
        return hash((self.i, self.j))

def solve():
    """
    Finds the number of power subgroups in the generalized quaternion group of size 128.
    """
    print("This script calculates the number of power subgroups in Q128.")
    
    # Generate all 128 elements of the group G = Q128
    group_G = []
    for i in range(64):
        group_G.append(Q128Element(i, 0))  # Elements of the form x^i
    for i in range(64):
        group_G.append(Q128Element(i, 1))  # Elements of the form x^i*y

    # The exponent of Q128 is 64. We check n from 0 to 64.
    exponent = 64
    power_subgroups = set()
    subgroup_generators = {}

    print(f"We compute the set G^n for n from 0 to {exponent} to find all distinct power subgroups.")
    
    for n in range(exponent + 1):
        g_n_set = set()
        for g in group_G:
            g_n_set.add(g**n)
        
        # Use frozenset as a key for dictionary and member of set
        f_g_n_set = frozenset(g_n_set)
        
        if f_g_n_set not in power_subgroups:
            subgroup_generators[f_g_n_set] = n
            power_subgroups.add(f_g_n_set)
    
    num_power_subgroups = len(power_subgroups)
    
    # Based on mathematical analysis, the subgroups are G, and the cyclic subgroups <x^(2^k)>
    descriptions_theory = [
        "The group G = Q128 itself",
        "The cyclic subgroup <x^2>",
        "The cyclic subgroup <x^4>",
        "The cyclic subgroup <x^8>",
        "The cyclic subgroup <x^16>",
        "The cyclic subgroup <x^32>",
        "The trivial subgroup {1}",
    ]

    print("\nFound the following distinct power subgroups, identified by the smallest non-negative integer n that generates them:")
    sorted_subgroups = sorted(list(power_subgroups), key=len, reverse=True)
    
    final_equation_parts = []
    for i, sg in enumerate(sorted_subgroups):
        n = subgroup_generators[sg]
        order = len(sg)
        final_equation_parts.append(str(order))
        desc = descriptions_theory[i] if i < len(descriptions_theory) else f"A subgroup"
        print(f"- G^{n:<2}: Subgroup of order {order:<3} ({desc})")

    print(f"\nThe different orders of the power subgroups are: {', '.join(sorted(final_equation_parts, key=int, reverse=True))}.")
    print(f"\nTotal number of power subgroups in Q128 is: {num_power_subgroups}")

if __name__ == '__main__':
    solve()