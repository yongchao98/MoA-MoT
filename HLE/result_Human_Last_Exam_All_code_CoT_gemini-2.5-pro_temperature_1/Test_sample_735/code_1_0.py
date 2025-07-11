import itertools

def is_product_free(n, s_tuple):
    """
    Checks if a set S (represented by s_tuple) is product-free in the
    cyclic group C_n (represented additively as Z_n).
    """
    s_set = set(s_tuple)
    for s1 in s_set:
        for s2 in s_set:
            # Using additive notation for the cyclic group C_n
            prod = (s1 + s2) % n
            if prod in s_set:
                return False
    return True

def has_maximal_product_free_set_of_size_3(n):
    """
    Checks if the cyclic group C_n has a maximal by inclusion
    product-free set of size 3.
    """
    group_elements = list(range(n))

    # Iterate through all subsets of size 3
    for s_tuple in itertools.combinations(group_elements, 3):
        if is_product_free(n, s_tuple):
            # The set s_tuple is product-free. Now check for maximality.
            # A set is maximal if adding any other element breaks the product-free property.
            is_maximal = True
            s_set = set(s_tuple)
            other_elements = [g for g in group_elements if g not in s_set]
            
            for g in other_elements:
                # Create a new set s_prime containing s_tuple and g
                s_prime_tuple = s_tuple + (g,)
                if is_product_free(n, s_prime_tuple):
                    # We found a larger product-free set, so s_tuple is not maximal.
                    is_maximal = False
                    break # No need to check other elements for this s_tuple
            
            if is_maximal:
                # We found a maximal product-free set of size 3.
                # We can stop searching and confirm this group qualifies.
                print(f"C_{n}: Found a maximal product-free set of size 3: {s_tuple}")
                return True

    print(f"C_{n}: No maximal product-free set of size 3 was found.")
    return False

def main():
    """
    Main function to check the list of cyclic groups from the literature.
    """
    # List of cyclic groups from the Giudici and Hart (2020) paper.
    cyclic_groups_to_check = [5, 7, 8, 9, 10, 11, 12, 13, 14]
    
    print("Checking for maximal product-free sets of size 3 in cyclic groups...\n")
    
    valid_groups = []
    for n in cyclic_groups_to_check:
        if has_maximal_product_free_set_of_size_3(n):
            valid_groups.append(f"C_{n}")

    print("\n--- Summary ---")
    print(f"The following cyclic groups from the list were confirmed to have a maximal product-free set of size 3:")
    print(", ".join(valid_groups))


if __name__ == "__main__":
    main()