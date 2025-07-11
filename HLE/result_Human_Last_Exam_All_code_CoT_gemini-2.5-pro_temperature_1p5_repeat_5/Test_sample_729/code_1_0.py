import math

class Q128Element:
    """
    Represents an element of the generalized quaternion group Q_128.
    An element is of the form x^i * y^j, where n=32 for Q_128.
    x^(2n)=1, y^2=x^n, y*x*y^(-1) = x^(-1)
    Here, 2n = 64, n = 32.
    """
    def __init__(self, i, j):
        # x^i y^j
        self.i = i % 64
        self.j = j % 2
        self.n = 32

    def __mul__(self, other):
        # (x^i1 y^j1) * (x^i2 y^j2)
        if self.j == 0:
            # x^i1 * x^i2 * y^j2 = x^(i1+i2) * y^j2
            return Q128Element(self.i + other.i, other.j)
        else: # self.j == 1
            # x^i1 * y * x^i2 * y^j2
            # using y * x^k = x^(-k) * y
            # x^i1 * x^(-i2) * y * y^j2 = x^(i1-i2) * y^(1+j2)
            new_i = self.i - other.i
            new_j = 1 + other.j
            if new_j == 2:
                # y^2 = x^n
                new_i += self.n
                new_j = 0
            return Q128Element(new_i, new_j)

    def __pow__(self, k):
        # Fast exponentiation
        if k == 0:
            return Q128Element(0, 0) # Identity element
        if k < 0:
            raise ValueError("Negative power not implemented for this group.")
            
        res = Q128Element(0, 0) # Start with identity
        base = self
        while k > 0:
            if k % 2 == 1:
                res = res * base
            base = base * base
            k //= 2
        return res
        
    def __eq__(self, other):
        return self.i == other.i and self.j == other.j

    def __hash__(self):
        return hash((self.i, self.j))

    def __repr__(self):
        if self.j == 0:
            return f"x^{self.i}"
        else:
            if self.i == 0:
                return "y"
            return f"x^{self.i}y"

def get_power_subgroups_count():
    """
    Calculates the number of power subgroups for Q128.
    """
    # Generate all elements of Q128
    G = []
    for i in range(64):
        for j in range(2):
            G.append(Q128Element(i, j))

    # The exponent of Q128 is 64. We only need to check divisors of 64.
    divisors_of_exponent = [1, 2, 4, 8, 16, 32, 64]
    
    found_subgroups = set()

    print("Finding power subgroups for Q_128:")
    print("-" * 40)
    
    total = 0
    for k in divisors_of_exponent:
        # Calculate G^k = {g^k for g in G}
        G_k = frozenset({g ** k for g in G})
        
        # Add the frozenset of elements to our set of found subgroups.
        # This automatically handles duplicates.
        if G_k not in found_subgroups:
            total += 1
            print(f"For k = {k}, a new power subgroup of size {len(G_k)} was found.")
            found_subgroups.add(G_k)
        else:
            print(f"For k = {k}, the power subgroup of size {len(G_k)} has already been found.")

    print("-" * 40)
    print(f"The total number of distinct power subgroups is: {len(found_subgroups)}")
    return len(found_subgroups)

if __name__ == '__main__':
    get_power_subgroups_count()