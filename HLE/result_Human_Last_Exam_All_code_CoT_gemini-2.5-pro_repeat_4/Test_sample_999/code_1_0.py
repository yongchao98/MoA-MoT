import numpy as np

def quantize_and_dequantize_demo(bits):
    """
    Demonstrates the process of uniform affine quantization and dequantization.
    This is the basic principle behind formats like INT8 and INT4.
    """
    print(f"--- Demonstrating {bits}-bit Uniform Quantization ---")

    # 1. Define the quantization integer range based on the number of bits.
    # For signed integers, the range is [-2^(n-1), 2^(n-1) - 1].
    if bits == 8:
        qmin, qmax = -128, 127
        dtype = np.int8
    elif bits == 4:
        qmin, qmax = -8, 7
        # NumPy doesn't have a native int4 type, so we use int8 to store the values
        # but ensure they are clamped within the 4-bit range.
        dtype = np.int8
    else:
        raise ValueError("Only 8 and 4 bits are supported in this demo.")

    # 2. Create a sample floating-point tensor (e.g., weights of a small layer).
    np.random.seed(42)
    float_tensor = (np.random.rand(2, 4) - 0.5) * 20  # Range approx [-10, 10]
    print("Original FP32 Tensor:\n", float_tensor)

    # 3. Calculate the scale and zero-point for quantization (per-tensor affine quantization).
    # Scale: maps the float range to the integer range.
    # Zero-point: ensures that the real value 0.0 maps to an integer.
    fmin, fmax = float_tensor.min(), float_tensor.max()
    
    # The scale is the change in float value for each integer step.
    scale = (fmax - fmin) / (qmax - qmin)
    
    # The zero point ensures the quantized range covers the original float range.
    # It's an integer.
    # Equation: zero_point = qmin - fmin / scale
    zero_point_float = qmin - fmin / scale
    zero_point = int(np.round(np.clip(zero_point_float, qmin, qmax)))

    print(f"\nCalculated Quantization Parameters (Equation Components):")
    print(f"Scale (S) = {scale:.4f}")
    print(f"Zero Point (Z) = {zero_point}")

    # 4. Quantize the tensor using the final equation: x_q = round(x / S) + Z
    quantized_tensor_float = np.round(float_tensor / scale) + zero_point
    quantized_tensor = np.clip(quantized_tensor_float, qmin, qmax).astype(dtype)
    print("\nQuantized INT Tensor (x_q):\n", quantized_tensor)

    # 5. Dequantize the tensor back to floating-point.
    # Final equation: x_approx = (x_q - Z) * S
    dequantized_tensor = (quantized_tensor.astype(np.float32) - zero_point) * scale
    print("\nDequantized FP32 Tensor (x_approx):\n", dequantized_tensor)

    # 6. Calculate and show the quantization error.
    error = np.mean((float_tensor - dequantized_tensor)**2)
    print(f"\nMean Squared Error between Original and Dequantized: {error:.4f}\n")


# Run the demo for INT8 as mentioned in the options
quantize_and_dequantize_demo(bits=8)

# Run the demo for INT4 as mentioned in the options
quantize_and_dequantize_demo(bits=4)