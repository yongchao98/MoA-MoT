import numpy as np

# Step 1: Define group structures and operations

# D = (C_2)^2 as vectors over F_2
D_elements = [(0, 0), (0, 1), (1, 0), (1, 1)]

def add_d(d1, d2):
    """Addition in D."""
    return ((d1[0] + d2[0]) % 2, (d1[1] + d2[1]) % 2)

# S = 3^{1+2}_+
# Elements are (a, b, c) representing x^a y^b z^c
S_elements = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]

def mult_s(s1, s2):
    """Multiplication in S."""
    a1, b1, c1 = s1
    a2, b2, c2 = s2
    # Using relation yx = xzy = x y z^{-1} -> y^b x^a = x^a y^b z^{-ab}
    a = (a1 + a2) % 3
    b = (b1 + b2) % 3
    c = (c1 + c2 - b1 * a2) % 3
    return (a, b, c)

def inv_s(s):
    """Inverse in S."""
    a, b, c = s
    inv_a = (-a) % 3
    inv_b = (-b) % 3
    inv_c = (-c + inv_b * a) % 3
    return (inv_a, inv_b, inv_c)

# G = D . S
G_elements = [(d, s) for d in D_elements for s in S_elements]

# Step 2: Implement the action phi: S -> Aut(D)
# Kernel K = <x, z> is the subgroup where b=0
K = {s for s in S_elements if s[1] == 0}

# Action is induced by S/K = C_3, map s=(a,b,c) -> b
# We need a non-trivial homomorphism C_3 -> Aut(D)
# Let's pick an order 3 matrix in GL(2, F_2)
A = np.array([[0, 1], [1, 1]])
A_sq = (A @ A) % 2
Id = np.identity(2, dtype=int)

# Cache powers of A
A_powers = {0: Id, 1: A, 2: A_sq}

def phi(s, d):
    """Action of s on d."""
    b = s[1]
    mat = A_powers[b]
    d_vec = np.array(d)
    res_vec = (mat @ d_vec) % 2
    return (res_vec[0], res_vec[1])

# Define G operations
def mult_g(g1, g2):
    """Multiplication in G."""
    d1, s1 = g1
    d2, s2 = g2
    s_part = mult_s(s1, s2)
    d_part = add_d(d1, phi(s1, d2))
    return (d_part, s_part)

def inv_g(g):
    """Inverse in G."""
    d, s = g
    s_inv = inv_s(s)
    # g^-1 = (s^-1, d^-1) in opposite group, here (phi(s^-1)(-d), s^-1)
    d_inv = phi(s_inv, d) # -d = d in char 2
    return (d_inv, s_inv)

# Step 3: Identify odd-order elements
identity_g = ((0, 0), (0, 0, 0))
odd_order_elements = set()
for g in G_elements:
    # An element has odd order iff g^3 = 1, since |s|=1 or 3 and |d|=1 or 2
    g2 = mult_g(g, g)
    g3 = mult_g(g2, g)
    s_part_g3 = g3[1]
    # s^3 is always 1 in S=3^{1+2}_+
    # We check the d component of g^3
    # g^3 = (d + phi(s)(d) + phi(s^2)(d), s^3)
    # The s part is always (0,0,0). So g3[1] is always identity.
    # We only need to check the d part.
    if g3[0] == (0, 0):
        odd_order_elements.add(g)

# Step 4: Partition the odd-order elements
odd_order_s_in_K = {g for g in odd_order_elements if g[1] in K}
odd_order_s_out_K = {g for g in odd_order_elements if g[1] not in K}

# Step 5: Function to count conjugacy classes
def count_conjugacy_classes(elements_subset, full_group):
    """Counts conjugacy classes of G within a given subset."""
    if not elements_subset:
        return 0
    
    remaining_elements = set(elements_subset)
    num_classes = 0
    
    while remaining_elements:
        num_classes += 1
        # Pick an element and compute its class
        g = remaining_elements.pop()
        current_class = {g}
        for h in full_group:
            h_inv = inv_g(h)
            conjugate = mult_g(mult_g(h, g), h_inv)
            current_class.add(conjugate)
        
        # Remove all elements of this class from the remaining set
        remaining_elements -= current_class
        
    return num_classes

# Step 6: Calculate and present the final result
num_classes_in_K = count_conjugacy_classes(odd_order_s_in_K, G_elements)
num_classes_out_K = count_conjugacy_classes(odd_order_s_out_K, G_elements)
total_blocks = num_classes_in_K + num_classes_out_K

print("The number of blocks is the number of 2'-conjugacy classes.")
print("This is calculated by counting the conjugacy classes of elements of odd order.")
print("We partition these elements based on whether their S-component is in K.")
print(f"Number of classes with s-component in K: {num_classes_in_K}")
print(f"Number of classes with s-component not in K: {num_classes_out_K}")
print(f"Total number of blocks = {num_classes_in_K} + {num_classes_out_K} = {total_blocks}")

print("\n<<<11>>>")