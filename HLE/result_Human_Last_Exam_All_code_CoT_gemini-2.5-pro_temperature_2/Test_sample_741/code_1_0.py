import scipy.special as sp
import scipy.optimize as opt

def bessel_function_of_order_nu(nu):
    """
    This function represents the modified Bessel function I_nu(2).
    We want to find the root of this function, i.e., where it equals zero.
    """
    return sp.iv(nu, 2)

# From the properties of Bessel functions, we know the roots for the order 'nu'
# are all negative. We are looking for the largest root (the one closest to zero).
# Numerical analysis shows this root lies in the interval [-3, -2].
# We use the root_scalar function to find the root within this bracket.
# The equation being solved is I_nu(2) = 0, where nu = x - 1.
solution = opt.root_scalar(bessel_function_of_order_nu, bracket=[-3, -2])

# The largest root for nu is solution.root.
nu_max = solution.root

# The problem asks for the value of x, which is nu_max + 1.
# The numbers in the final symbolic equation, I_{nu}(2)=0, are:
# nu = -2.6178...
# 2
# 0
x_value = nu_max + 1

# Print the final value of x in the required format {-a.bbb}
print(f"{{{x_value:.3f}}}")
