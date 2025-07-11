import collections
from itertools import combinations

class FiniteGroup:
    """A class to represent a finite group using its Cayley table."""
    def __init__(self, name, elements, cayley_table):
        self.name = name
        self.elements = frozenset(elements)
        self.cayley_table = cayley_table
        # Find the identity element
        for e in self.elements:
            if all(self.product(e, g) == g for g in self.elements):
                self.identity_elt = e
                break

    def product(self, a, b):
        """Calculates the product of two elements in the group."""
        return self.cayley_table[a][b]

    def __repr__(self):
        return self.name

def is_product_free(s, group):
    """Checks if a set s is product-free within a group."""
    for x in s:
        for y in s:
            if group.product(x, y) in s:
                return False
    return True

def has_maximal_product_free_set_of_size_2(group):
    """
    Checks if a group has a maximal by inclusion product-free set of size 2.
    """
    potential_elements = list(group.elements - {group.identity_elt})
    if len(potential_elements) < 2:
        return False

    # Iterate through all subsets of size 2
    for s_tuple in combinations(potential_elements, 2):
        s = set(s_tuple)

        # 1. Check if S is product-free
        if not is_product_free(s, group):
            continue

        # 2. If product-free, check for maximality
        is_maximal = True
        other_elements = group.elements - s
        for g in other_elements:
            t = s.union({g})
            if is_product_free(t, group):
                # Found a larger product-free set, so s is not maximal
                is_maximal = False
                break
        
        if is_maximal:
            # This group has the property
            return True
            
    return False

# --- Group Definitions ---

def make_zn(n):
    """Creates the cyclic group Z_n."""
    elements = list(range(n))
    table = {i: {j: (i + j) % n for j in elements} for i in elements}
    return FiniteGroup(f"Z_{n}", elements, table)

def make_v4():
    """Creates the Klein four-group V_4."""
    elements = ['e', 'a', 'b', 'c']
    table = {
        'e': {'e':'e', 'a':'a', 'b':'b', 'c':'c'},
        'a': {'e':'a', 'a':'e', 'b':'c', 'c':'b'},
        'b': {'e':'b', 'a':'c', 'b':'e', 'c':'a'},
        'c': {'e':'c', 'a':'b', 'b':'a', 'c':'e'},
    }
    return FiniteGroup("V_4", elements, table)

def make_d2n(n):
    """Creates the dihedral group D_2n of order 2n."""
    elements = [f'r{i}' for i in range(n)] + [f's{i}' for i in range(n)]
    table = {el1: {el2: "" for el2 in elements} for el1 in elements}
    for i in range(n):
        for j in range(n):
            table[f'r{i}'][f'r{j}'] = f'r{(i+j)%n}'      # r^i * r^j = r^(i+j)
            table[f'r{i}'][f's{j}'] = f's{(i+j)%n}'      # r^i * sr^j = s r^(-i)r^j = sr^(j-i) is wrong. s_j = r^j s. r^i s_j = r^(i+j)s
            table[f's{i}'][f'r{j}'] = f's{(i-j+n)%n}'    # s_i r^j = r^i s r^j = r^i r^(-j) s = r^(i-j) s
            table[f's{i}'][f's{j}'] = f'r{(i-j+n)%n}'    # s_i s_j = r^i s r^j s = r^i r^(-j) s s = r^(i-j)
    
    # Use more standard names like e, r, s, sr, etc.
    old_to_new = {f'r{i}':(f'r{i}' if i > 0 else 'e') for i in range(n)}
    old_to_new.update({f's{i}':(f's' if i==0 else f'sr{i}') for i in range(n)})
    
    new_elements = sorted(list(old_to_new.values()), key=lambda x: (len(x), x))
    new_table = {old_to_new[e1]: {old_to_new[e2]: old_to_new[table[e1][e2]] for e2 in elements} for e1 in elements}
    return FiniteGroup(f"D_{2*n}", new_elements, new_table)

if __name__ == '__main__':
    # List of non-isomorphic groups to check
    groups_to_check = []
    # Cyclic groups
    for n in range(3, 11):
        groups_to_check.append(make_zn(n))
    # Other common small groups
    groups_to_check.append(make_v4())
    groups_to_check.append(make_d2n(3)) # S_3 is isomorphic to D_6
    groups_to_check.append(make_d2n(4)) # D_8
    groups_to_check.append(make_d2n(5)) # D_10

    found_groups = []
    for group in groups_to_check:
        if has_maximal_product_free_set_of_size_2(group):
            found_groups.append(group.name)

    print(f"Found {len(found_groups)} groups with the specified property.")
    print("The groups are:", ", ".join(sorted(found_groups)))
    print("\nFinal Answer:")
    print(len(found_groups))