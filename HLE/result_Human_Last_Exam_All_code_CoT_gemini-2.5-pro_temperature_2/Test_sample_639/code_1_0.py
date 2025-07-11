import math

def find_packing_distances():
    """
    Calculates and prints the set of possible normalized distances 'r'
    between hard spheres packed on a 2D plane for r <= 3.

    The calculation is based on a triangular lattice model, which represents the
    densest packing of spheres in 2D. The normalized distance 'r' from a
    central sphere to any other sphere can be found using the formula:
    r^2 = n^2 + n*m + m^2, where n and m are integers defining the lattice position.

    This script finds all unique distances for r <= 3.
    """
    r_max = 3
    r_max_sq = r_max**2
    
    # Using a set to automatically store unique values of r^2
    distances_sq_set = set()

    # A search range for n and m from -3 to 3 is sufficient, as larger values
    # would result in r^2 > 9.
    search_limit = r_max
    
    for n in range(-search_limit, search_limit + 1):
        for m in range(-search_limit, search_limit + 1):
            # (n=0, m=0) is the reference sphere itself, so we skip it.
            if n == 0 and m == 0:
                continue

            r_sq = float(n**2 + n * m + m**2)

            if r_sq <= r_max_sq:
                distances_sq_set.add(r_sq)

    # Calculate r from the unique r^2 values and sort them.
    sorted_distances = sorted([math.sqrt(r_sq) for r_sq in distances_sq_set])
    
    print("The calculated set of normalized distances (r) is:")
    # "Remember in the final code you still need to output each number in the final equation!"
    # Interpreting "the final equation" as the final resulting set of numbers.
    # We print each value, formatted to two decimal places.
    for r in sorted_distances:
        # To make it look like an "equation", we can show the source of r
        r_val_string = f"{r:.2f}"
        if r_val_string == "1.00":
             print(f"r = sqrt(1^2 + 1*0 + 0^2) = sqrt(1.00) = {r_val_string}")
        elif r_val_string == "1.73":
             print(f"r = sqrt(1^2 + 1*1 + 1^2) = sqrt(3.00) = {r_val_string}")
        elif r_val_string == "2.00":
             print(f"r = sqrt(2^2 + 2*0 + 0^2) = sqrt(4.00) = {r_val_string}")
        elif r_val_string == "2.65":
             print(f"r = sqrt(2^2 + 2*1 + 1^2) = sqrt(7.00) = {r_val_string}")
        elif r_val_string == "3.00":
             print(f"r = sqrt(3^2 + 3*0 + 0^2) = sqrt(9.00) = {r_val_string}")

if __name__ == "__main__":
    find_packing_distances()