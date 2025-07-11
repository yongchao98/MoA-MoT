import sys
from codecs import encode

def formulate_it3_vertical_slice():
    """
    This function formulates and prints the mathematical equation for the
    vertical cross-section of a Gaussian-based Interval Type-3 Membership Function.
    """

    # For broader terminal compatibility, we'll use ASCII/Unicode representations.
    # The 'overline' character combines a letter with the U+0305 combining overline.
    mu_overline = 'μ' + u'\u0305'
    sigma_overline = 'σ' + u'\u0305'
    c_overline = 'c' + u'\u0305'

    # The equation describes the upper bound of the FOU for the vertical slice.
    # It is a function of the secondary variable 'u' for a fixed primary input 'x_p'.
    # The key feature of the IT3 MF is that the Gaussian parameters 'c' and 'σ'
    # are themselves functions of 'u'.
    
    print("The mathematical formulation characterizing the upper bound of a vertical cross-section of a Gaussian IT3 MF is given by:")
    print("-" * 80)
    
    # Constructing and printing the final equation string.
    # We explicitly show each variable and number as requested.
    equation = (
        f"{mu_overline}_vc(u) = exp( -0.5 * ( (x_p - {c_overline}(u)) / {sigma_overline}(u) )^2 )"
    )
    print(equation)
    
    print("-" * 80)
    print("\nWhere:")
    print(f"  - {mu_overline}_vc(u): Represents the upper membership grade for the vertical cross-section. It is a function of the secondary variable 'u'.")
    print("  - x_p: A specific, fixed value of the primary input variable 'x'.")
    print("  - u: The secondary variable, representing the vertical axis of the secondary membership function, typically in the range [0, 1].")
    print("  - exp(): The standard exponential function.")
    print("  - The number -0.5 and the exponent 2 define the Gaussian form.")
    print(f"  - {c_overline}(u): The center (mean) of the Gaussian function. In an IT3 MF, this center is a function of 'u'.")
    print(f"  - {sigma_overline}(u): The spread (standard deviation) of the Gaussian function, which is also a function of 'u'.")
    print("\nThis equation defines the Upper Membership Function (UMF) of the Footprint of Uncertainty (FOU) for the Interval Type-2 set that results from the vertical slice.")
    print("A similar equation with different parameter functions would define the Lower Membership Function (LMF).")


if __name__ == "__main__":
    # Handle potential encoding issues in different terminals
    if sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
        sys.stdout.reconfigure(encoding='utf-8')
    formulate_it3_vertical_slice()

<<<
The mathematical formulation characterizing the upper bound of a vertical cross-section of a Gaussian IT3 MF is given by:
--------------------------------------------------------------------------------
μ̄_vc(u) = exp( -0.5 * ( (x_p - c̄(u)) / σ̄(u) )^2 )
--------------------------------------------------------------------------------

Where:
  - μ̄_vc(u): Represents the upper membership grade for the vertical cross-section. It is a function of the secondary variable 'u'.
  - x_p: A specific, fixed value of the primary input variable 'x'.
  - u: The secondary variable, representing the vertical axis of the secondary membership function, typically in the range [0, 1].
  - exp(): The standard exponential function.
  - The number -0.5 and the exponent 2 define the Gaussian form.
  - c̄(u): The center (mean) of the Gaussian function. In an IT3 MF, this center is a function of 'u'.
  - σ̄(u): The spread (standard deviation) of the Gaussian function, which is also a function of 'u'.

This equation defines the Upper Membership Function (UMF) of the Footprint of Uncertainty (FOU) for the Interval Type-2 set that results from the vertical slice.
A similar equation with different parameter functions would define the Lower Membership Function (LMF).
>>>