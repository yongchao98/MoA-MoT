import numpy as np

# A color is unrepresentable in sRGB if, after conversion, any of its R, G, B
# components are outside the range [0.0, 1.0]. This script performs the standard
# conversion from ROMM RGB to sRGB and checks this condition for each color.

# 1. Define constants and matrices for the color space conversion.
# The ROMM RGB color space (also known as ProPhoto RGB) uses a gamma of 1.8.
ROMM_GAMMA = 1.8

# The full transformation from linear ROMM RGB to linear sRGB is a combination
# of three matrix operations.

# Matrix a: Convert from ProPhoto RGB (D50 white point) to the CIE XYZ color space.
M_prophoto_to_xyz = np.array([
    [0.7976749, 0.1351917, 0.0313534],
    [0.2880402, 0.7118741, 0.0000857],
    [0.0000000, 0.0000000, 0.8252100]
])

# Matrix b: Apply Bradford chromatic adaptation from D50 (ProPhoto's white point)
# to D65 (sRGB's white point).
M_cat_d50_to_d65 = np.array([
    [ 0.9555766, -0.0230393, 0.0631636],
    [-0.0282895,  1.0099416, 0.0210077],
    [ 0.0122982, -0.0204830, 1.3299098]
])

# Matrix c: Convert from the CIE XYZ color space to linear sRGB.
M_xyz_to_srgb = np.array([
    [ 3.2404542, -1.5371385, -0.4985314],
    [-0.9692660,  1.8760108,  0.0415560],
    [ 0.0556434, -0.2040259,  1.0572252]
])

# We multiply these matrices together to get a single, final transformation matrix.
M_ROMM_TO_SRGB = M_xyz_to_srgb @ M_cat_d50_to_d65 @ M_prophoto_to_xyz

# 2. Define the input ROMM RGB colors to be tested.
colors = {
    1: {"name": "RGB(0, 0, 1)", "values": np.array([0.0, 0.0, 1.0])},
    2: {"name": "RGB(0, 1, 0)", "values": np.array([0.0, 1.0, 0.0])},
    3: {"name": "RGB(0, 0.5, 0.6)", "values": np.array([0.0, 0.5, 0.6])},
    4: {"name": "RGB(0.4, 0.5, 0.6)", "values": np.array([0.4, 0.5, 0.6])},
    5: {"name": "RGB(1, 1, 1)", "values": np.array([1.0, 1.0, 1.0])}
}

out_of_gamut_indices = []

print("--- Checking which ROMM RGB colors can be represented in sRGB ---")

# 3. Process each color.
for index, color_data in colors.items():
    romm_rgb = color_data["values"]
    
    # Step A: Linearize the ROMM RGB values by applying the inverse gamma.
    linear_romm_rgb = np.power(romm_rgb, ROMM_GAMMA)
    
    # Step B: Convert from linear ROMM to linear sRGB using the final matrix.
    linear_srgb = M_ROMM_TO_SRGB @ linear_romm_rgb
    
    # Step C: Check if any component is outside the valid [0, 1] range.
    # A small tolerance (-1e-6) is used to account for floating-point inaccuracies.
    is_out_of_gamut = np.any(linear_srgb < -1e-6) or np.any(linear_srgb > 1 + 1e-6)
    
    # Output the details of the calculation.
    print(f"\nColor {index}: {color_data['name']}")
    print(f"  - Calculated linear sRGB values: R={linear_srgb[0]:.4f}, G={linear_srgb[1]:.4f}, B={linear_srgb[2]:.4f}")

    if is_out_of_gamut:
        out_of_gamut_indices.append(index)
        print("  - Result: Cannot be represented in sRGB (values are outside the [0, 1] range).")
    else:
        print("  - Result: Can be represented in sRGB.")
        
# 4. Phrase and print the final answer.
if not out_of_gamut_indices:
    result_str = "none cannot"
else:
    # Sort the indices for the final answer.
    result_str = ", ".join(map(str, sorted(out_of_gamut_indices)))

print("\n--- Final Answer ---")
print(f"The numbers of the colors that cannot be represented by an sRGB hex code are: {result_str}")
<<<1, 2, 3, 5>>>