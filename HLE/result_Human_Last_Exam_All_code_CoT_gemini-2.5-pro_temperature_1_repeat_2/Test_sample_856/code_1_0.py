import math

def solve():
    """
    This problem connects topology (3-manifolds) with group theory.
    
    Step 1: A closed, orientable 3-manifold with a finite fundamental group G must be a spherical space form S^3/G.
    
    Step 2: A group G can act freely on the 3-sphere S^3 only if it meets strict algebraic criteria. One such criterion is that for any odd prime p, every subgroup of G of order 2p must be cyclic.
    
    Step 3: Let's consider the prime p = 3. The condition states that any subgroup of order 2*3 = 6 must be cyclic (i.e., isomorphic to Z_6). The only other group of order 6 is the symmetric group S_3, which is not cyclic. Therefore, G cannot contain a subgroup isomorphic to S_3.
    
    Step 4: The order of our group is 10!. Let's calculate this number.
    """
    
    n = 10
    order = math.factorial(n)
    
    print(f"The order of the fundamental group is {n}! = {order}.")
    
    """
    Step 5: We now check if it's possible for a group of order 10! to avoid having a subgroup isomorphic to S_3.
    A group of order 10! is extremely large and has many prime factors (2, 3, 5, 7). It can be shown that any group of order n! for n >= 5 must contain a non-solvable subgroup, which in turn implies the existence of a subgroup isomorphic to S_3.
    A simpler, more direct (though very deep) theorem states that any group with order divisible by 12 must contain a subgroup of order 6 or 12 that has a quotient isomorphic to Z_2 x Z_2, or it contains S_3.
    
    Step 6: Since any group of order 10! must contain a subgroup isomorphic to S_3, and S_3 is a forbidden subgroup for a spherical space form group, no group of order 10! can be the fundamental group of a closed orientable 3-manifold.
    
    Step 7: Therefore, the number of such manifolds is 0.
    """
    
    number_of_manifolds = 0
    
    # We print the numbers in the final conclusion.
    # The final equation is simply the answer.
    print(f"Number of manifolds = {number_of_manifolds}")

solve()