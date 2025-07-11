import numpy as np

def demonstrate_quantization_scheme():
    """
    Demonstrates the per-token/per-channel quantization scheme for a linear layer.
    This illustrates the concept behind statement B.
    """
    # Define matrix dimensions
    T = 8      # Number of tokens (rows in activation matrix)
    C_i = 64   # Input features / channels
    C_o = 128  # Output features / channels

    # --- Step 1: Create random FP32 data ---
    np.random.seed(42)
    X_fp32 = (np.random.randn(T, C_i) * 5).astype(np.float32)
    W_fp32 = (np.random.randn(C_i, C_o) * 0.1).astype(np.float32)

    # --- Step 2: Perform ground truth calculation in FP32 ---
    Y_fp32 = X_fp32 @ W_fp32

    print("This code demonstrates the quantization process described in statement B.")
    print(f"Goal: Compute Y = XW, where X is ({T}x{C_i}) and W is ({C_i}x{C_o}).")
    print("-" * 60)

    # --- Step 3: Quantize tensors to INT8 ---
    # Per-token (per-row) symmetric quantization for activations X
    # scale = max(abs(row)) / 127
    s_X = np.max(np.abs(X_fp32), axis=1, keepdims=True) / 127.0
    X_int8 = np.round(np.clip(X_fp32 / (s_X + 1e-9), -127, 127)).astype(np.int8)

    # Per-channel (per-column) symmetric quantization for weights W
    # scale = max(abs(column)) / 127
    s_W = np.max(np.abs(W_fp32), axis=0, keepdims=True) / 127.0
    W_int8 = np.round(np.clip(W_fp32 / (s_W + 1e-9), -127, 127)).astype(np.int8)

    print("Step A: Quantization (FP32 -> INT8)")
    print(f"Activation scales 's_X' shape (one per token): {s_X.shape}")
    print(f"Weight scales 's_W' shape (one per output channel): {s_W.shape}")
    print("-" * 60)

    # --- Step 4: Perform matrix multiplication with INT8 tensors ---
    # On a GPU, this is a highly accelerated INT8 GEMM. Result is INT32.
    Y_int_gemm = X_int8.astype(np.int32) @ W_int8.astype(np.int32)
    
    print("Step B: INT8 Matrix Multiplication (GEMM)")
    print(f"Resulting integer matrix 'Y_int_gemm' shape: {Y_int_gemm.shape}")
    print("-" * 60)

    # --- Step 5: De-quantize the result ---
    # The combined scale is the outer product of the two scale vectors
    # Broadcasting: (T, 1) * (1, Co) -> (T, Co)
    combined_scale = s_X * s_W
    Y_quantized_fp32 = Y_int_gemm.astype(np.float32) * combined_scale

    print("Step C: De-quantization (INT32 -> FP32)")
    print("The de-quantization equation is: Y_quantized[i,j] = Y_int_gemm[i,j] * s_X[i] * s_W[j]")
    
    # Show the calculation for the first element Y[0,0] as an example
    y_int_00 = Y_int_gemm[0, 0]
    s_x_0 = s_X[0, 0]
    s_w_0 = s_W[0, 0]
    y_quant_00 = Y_quantized_fp32[0, 0]
    
    print("\nExample calculation for the element at [0, 0]:")
    print(f"  Y_quantized[0,0] = {y_int_00} * {s_x_0:.6f} * {s_w_0:.6f}")
    print(f"  Result: {y_quant_00:.6f}")
    print("-" * 60)

    # --- Step 6: Compare results ---
    mse = np.mean((Y_fp32 - Y_quantized_fp32)**2)
    print("Step D: Verification")
    print("Comparing the first 4x4 block of the matrices:")
    print("Original FP32 result:\n", Y_fp32[:4, :4])
    print("\nReconstructed result from INT8 GEMM:\n", Y_quantized_fp32[:4, :4])
    print(f"\nMean Squared Error between original and quantized results: {mse:.8f}")
    print("The low error shows the method is effective.")

if __name__ == '__main__':
    demonstrate_quantization_scheme()