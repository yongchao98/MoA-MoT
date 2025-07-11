def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def find_filled_groups(q_max, m_max):
    """
    Finds and describes nonabelian filled groups of order 2*q^m for
    odd primes q up to q_max and natural numbers m up to m_max.
    """
    print(f"Searching for nonabelian filled groups of order 2*q^m for q <= {q_max}, m <= {m_max}")
    
    found_groups = []
    
    # Iterate through possible values of q (odd primes)
    for q in range(3, q_max + 1):
        if not is_prime(q) or q == 2:
            continue
            
        # Iterate through possible values of m (natural numbers)
        for m in range(1, m_max + 1):
            
            # The condition derived is that m must be a positive even integer.
            if m > 0 and m % 2 == 0:
                group_order_val = 2 * (q**m)
                order_str = f"2 * {q}^{m} = {group_order_val}"
                description = (
                    f"For q={q}, m={m}: A nonabelian filled group of order {order_str} exists.\n"
                    f"It is the generalized dihedral group Dih((Z_{q})^{m}), which is the semidirect product\n"
                    f"(Z_{q})^{m} x| Z_2, where Z_q is the cyclic group of order {q} and the element of\n"
                    f"order 2 acts on (Z_{q})^{m} by sending each element to its inverse."
                )
                print("-" * 20)
                print(description)
                found_groups.append({'q': q, 'm': m})
            else:
                # Uncomment the following lines to see the cases that are not filled groups.
                # print("-" * 20)
                # print(f"For q={q}, m={m}: The nonabelian group of order 2*{q}^{m} is not a filled group according to the classification (m must be even).")
                pass
                
    if not found_groups:
        print("No nonabelian filled groups found for the given parameters.")

if __name__ == '__main__':
    # Example: Check for odd primes up to 7 and m up to 4.
    find_filled_groups(q_max=7, m_max=4)