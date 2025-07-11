import math

def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def demonstrate_construction():
    """
    Demonstrates the construction of an infinite set 'x' that is
    almost disjoint from a countable family of infinite sets 'S'.
    """
    # Let S = {s_0, s_1, s_2, ...} be a countable family of infinite sets.
    # We represent these infinite sets as boolean functions.
    # s_0: set of even numbers
    # s_1: set of multiples of 3
    # s_2: set of prime numbers
    # s_3: set of perfect squares
    # s_4: set of numbers ending in 7
    s = [
        lambda n: n % 2 == 0,
        lambda n: n % 3 == 0,
        lambda n: is_prime(n),
        lambda n: int(math.sqrt(n))**2 == n if n >= 0 else False,
        lambda n: n % 10 == 7
    ]
    
    num_sets_to_consider = len(s)
    num_elements_to_find = 15
    
    print(f"Demonstrating the construction for {num_elements_to_find} elements of x.")
    print("The family S consists of the following sets:")
    print("s_0: Even numbers")
    print("s_1: Multiples of 3")
    print("s_2: Prime numbers")
    print("s_3: Perfect squares")
    print("s_4: Numbers ending in 7")
    print("-" * 30)
    
    x = []
    x_previous = -1
    
    # Construct the set x = {x_0, x_1, x_2, ...}
    for k in range(num_elements_to_find):
        # We need to find x_k > x_{k-1} such that x_k is not in s_0, s_1, ..., s_k.
        # This means we must avoid the union of the first k+1 sets.
        
        candidate = x_previous + 1
        while True:
            is_in_union = False
            # Check against s_0, s_1, ..., s_k (or up to the number of sets defined)
            num_sets_to_check = min(k + 1, num_sets_to_consider)
            
            for i in range(num_sets_to_check):
                if s[i](candidate):
                    is_in_union = True
                    break
            
            if not is_in_union:
                # Found our element for x_k
                x_k = candidate
                x.append(x_k)
                x_previous = x_k
                print(f"Found x_{k} = {x_k}")
                break
            
            candidate += 1
            
    print("\n" + "-" * 30)
    print(f"The first {num_elements_to_find} elements of the constructed set x are:")
    print(f"x = {x}")
    print("-" * 30)
    
    # Verify that for each s_i, the intersection x_intersect_s_i is finite.
    print("Verifying the intersections |x intersect s_i| are finite:")
    for i in range(num_sets_to_consider):
        # By construction, x_k is not in s_i for all k >= i.
        # So the intersection can only contain {x_0, ..., x_{i-1}}.
        intersection = {val for val in x if s[i](val)}
        
        # Elements of x with index >= i
        tail_of_x = x[i:]
        tail_intersection = {val for val in tail_of_x if s[i](val)}

        print(f"\nIntersection of x with s_{i}:")
        print(f"x_intersect_s_{i} = {intersection if intersection else '{}'}")
        print(f"Size = {len(intersection)}")
        
        # The crucial check:
        print(f"Note: By construction, no element x_k for k >= {i} is in s_{i}.")
        print(f"Intersection with tail x_{'>='}{i} is: {tail_intersection if tail_intersection else '{}'}")
        

demonstrate_construction()