import sys

def is_prime(n):
    """A simple primality test function."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def analyze_filled_groups(q, m):
    """
    Identifies the nonabelian filled groups of order 2*q**m based on
    the mathematical characterization.
    
    Args:
        q: An odd prime number.
        m: A natural number.
    """
    # Step 1: Validate the input parameters
    if not (isinstance(q, int) and q > 2 and q % 2 != 0 and is_prime(q)):
        print(f"Error: The input q must be an odd prime number. Received q={q}.")
        return
    if not (isinstance(m, int) and m >= 1):
        print(f"Error: The input m must be a natural number (m >= 1). Received m={m}.")
        return

    # Step 2: Apply the derived condition (m must be 1)
    print(f"Analyzing nonabelian groups of order 2 * {q}^{m} = {2 * q**m}...")
    if m == 1:
        order = 2 * q
        print(f"For q={q} and m={m}, a nonabelian filled group exists.")
        print(f"The group is the Dihedral group of order {order}.")
        print("Final Equation: D_(2 * q)")
        print(f"For the given q, this is: D_(2 * {q}) = D_{order}")
    else:
        print(f"For q={q} and m={m}, no nonabelian filled groups of order {2 * q**m} exist.")
        print("This is because the condition for such a group to be 'filled' requires m=1.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <q> <m>")
        print("  q: an odd prime number (e.g., 3, 5, 7, ...)")
        print("  m: a natural number (e.g., 1, 2, 3, ...)")
        sys.exit(1)
    
    try:
        q_val = int(sys.argv[1])
        m_val = int(sys.argv[2])
        analyze_filled_groups(q_val, m_val)
    except ValueError:
        print("Error: q and m must be integers.")
        sys.exit(1)