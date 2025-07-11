import collections
from itertools import product

class AbelianGroup:
    """A class to represent a finite Abelian group as a product of cyclic groups."""
    def __init__(self, orders):
        self.orders = tuple(orders)
        self.size = 1
        for o in orders:
            self.size *= o
        self.dim = len(orders)
        self.identity = tuple([0] * self.dim)
        self._elements = list(self._generate_elements())

    def _generate_elements(self):
        """Generate all elements of the group."""
        ranges = [range(o) for o in self.orders]
        return list(product(*ranges))

    def elements(self):
        return self._elements

    def add(self, a, b):
        """Add two elements in the group."""
        return tuple([(a[i] + b[i]) % self.orders[i] for i in range(self.dim)])

    def double(self, a):
        """Double an element."""
        return tuple([(2 * a[i]) % self.orders[i] for i in range(self.dim)])

def get_partitions(n):
    """Get all integer partitions of n."""
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def get_abelian_group_structures(n):
    """Generate structures of all abelian groups of order n."""
    if n == 1:
        yield [1]
        return
    
    # Prime factorization of n
    factors = collections.defaultdict(int)
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    
    # Get partitions for each prime power factor
    prime_power_partitions = []
    for p, exponent in factors.items():
        partitions_for_p = []
        for part in get_partitions(exponent):
            partitions_for_p.append([p**i for i in part])
        prime_power_partitions.append(partitions_for_p)

    # Combine partitions using Cartesian product
    for combo in product(*prime_power_partitions):
        structure = []
        for part in combo:
            structure.extend(part)
        yield sorted(structure)

def find_solution():
    """Main function to find the smallest group size."""
    
    order_n = 1
    while True:
        order_n += 1
        if order_n % 2 != 0 and order_n > 1: # Skip odd orders
             continue
        if order_n in [2]: continue # Skip known failing cases

        for struct in get_abelian_group_structures(order_n):
            group = AbelianGroup(struct)
            
            # Find G_2 and 2G
            g2_set = {g for g in group.elements() if group.double(g) == group.identity}
            two_g_set = {group.double(g) for g in group.elements()}
            
            elements = [g for g in group.elements() if g != group.identity]
            
            # Backtracking search for maximal sum-free sets
            memo_max = set()
            
            # stack stores (current_set, potential_additions)
            stack = [(frozenset(), tuple(elements))]
            
            while stack:
                current_s, potential_adds = stack.pop()
                
                is_maximal = True
                for i in range(len(potential_adds)):
                    g = potential_adds[i]
                    
                    # Test if adding g preserves sum-freeness
                    next_s = current_s | {g}
                    is_sf = True
                    for x in current_s:
                        if group.add(x, g) in next_s:
                            is_sf = False
                            break
                    if is_sf and group.add(g, g) in next_s:
                        is_sf = False

                    if is_sf:
                        is_maximal = False
                        new_potential_adds = potential_adds[:i] + potential_adds[i+1:]
                        stack.append((next_s, new_potential_adds))

                if is_maximal and current_s:
                    if current_s in memo_max:
                        continue
                    memo_max.add(current_s)
                    
                    s = current_s
                    s_cap_2g = s.intersection(two_g_set)
                    k_s_size = len(s_cap_2g) * len(g2_set)

                    if k_s_size > 2 * len(s):
                        print(f"Found a solution in group G of order {group.size}")
                        print(f"Group Structure (Z_orders): {group.orders}")
                        print(f"Maximal sum-free set S: {s}")
                        print(f"|S| = {len(s)}")
                        print(f"k(S) = {{g in G | 2g in S}}")
                        print(f"|k(S)| = {k_s_size}")
                        print(f"Check: {k_s_size} > 2 * {len(s)} -> {k_s_size > 2 * len(s)}")
                        return group.size

find_solution()