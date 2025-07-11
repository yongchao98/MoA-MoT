import itertools

def get_limit(poset_objects, sets, maps):
    """
    Calculates the limit of a diagram F: J^op -> Set for a finite poset J.

    Args:
        poset_objects: A list of objects in the poset J.
        sets: A dictionary mapping each object j in J to a set F(j).
        maps: A dictionary mapping each morphism (i, j) where i <= j to a function F(j) -> F(i).
              The partial order is inferred from the keys of the maps dictionary.
    """
    limit_set = []
    
    # The limit is a subset of the Cartesian product of the sets.
    # We create the product space.
    set_list = [sets[i] for i in poset_objects]
    
    if any(not s for s in set_list):
        product_space = [] # Product is empty if any constituent set is empty.
    else:
        product_space = itertools.product(*set_list)

    for p_tuple in product_space:
        # p_tuple is an element of the product space, e.g., (x_0, x_1, ...).
        # We map it to a dictionary for easier lookup by object name.
        p_dict = {obj: val for obj, val in zip(poset_objects, p_tuple)}
        
        is_in_limit = True
        # Check the compatibility condition for all morphisms.
        # A morphism exists from j to i if (i, j) is a key in `maps`.
        for i, j in maps.keys():
            f_ij = maps[(i, j)]
            # The condition is f_ij(x_j) = x_i
            if f_ij(p_dict[j]) != p_dict[i]:
                is_in_limit = False
                break
        
        if is_in_limit:
            limit_set.append(p_tuple)
            
    return limit_set

# --- Main part of the script ---

# 1. Define a simple finite directed poset J = ({0, 1, 2, 3}, <=)
J_objects = sorted([0, 1, 2, 3])

# 2. Define a functor F that maps every object to a non-empty set.
#    To find the minimum size, we choose the smallest possible non-empty set: a singleton.
F_sets = {i: {0} for i in J_objects}

# 3. Define the surjective maps for the functor.
#    For any i <= j, we need a surjective map F(j) -> F(i).
#    Since F(j) = {0} and F(i) = {0}, the only possible map is 0 -> 0. This is surjective.
def constant_map_to_zero(x):
    return 0

F_maps = {}
for i in J_objects:
    for j in J_objects:
        if i <= j:
            F_maps[(i, j)] = constant_map_to_zero

# 4. Compute the limit for this example system.
limit = get_limit(J_objects, F_sets, F_maps)
limit_size = len(limit)

# 5. Print the results and the verification of the equations.
print("--- Analysis of the minimum limit size ---")
print("Based on the non-empty limit theorem, the size must be >= 1.")
print("We test if a size of 1 is possible with a concrete example.\n")

print("--- Example Calculation ---")
print(f"Poset J = {J_objects} with the usual order <=.")
print(f"Sets F(j) = {F_sets[0]} for all j in J.")
print("Maps f_ij: F(j) -> F(i) for i <= j are all the function `x -> 0`.\n")

print(f"The calculated limit set is: {limit}")
print(f"The size of the limit set is: {limit_size}\n")

print("--- Verification of the Limit Element(s) ---")
if not limit:
    print("The limit is empty, which contradicts the theory (and shouldn't happen here).")
else:
    for p_tuple in limit:
        print(f"Checking element {p_tuple}:")
        p_dict = {obj: val for obj, val in zip(J_objects, p_tuple)}
        for i, j in sorted(list(F_maps.keys())):
            map_ij = F_maps[(i, j)]
            x_j = p_dict[j]
            x_i_expected = p_dict[i]
            x_i_actual = map_ij(x_j)
            
            # This part prints each number in the final equation for each condition
            print(f"  Condition for i={i}, j={j}: f_ij(x_{j}) = x_{i}")
            print(f"    f_{i}{j}({x_j}) = {x_i_actual}  (Expected: {x_i_expected}) -> Match: {x_i_actual == x_i_expected}")
        print(f"Conclusion: Element {p_tuple} satisfies all conditions.\n")

print("--- Final Answer ---")
print("Since the limit size must be >= 1, and we constructed a valid system where the size is 1,")
print("the smallest possible size of the limit is 1.")
