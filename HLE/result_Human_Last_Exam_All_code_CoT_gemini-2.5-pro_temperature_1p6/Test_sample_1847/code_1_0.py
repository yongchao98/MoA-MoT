import math

def print_analytical_solution():
  """
  Prints the analytical solution for the density profile of the specified shock wave.
  """
  
  # The density profile rho'(xi) is rho(x)/rho_0, where rho_0 is the ambient density.
  # The position xi is x/L, where L is the ambient conductive length scale.
  # The equation is derived from the Navier-Stokes equations for the specific case of
  # a Mach sqrt(3) shock in a monatomic ideal gas with a Prandtl number of 3/4.
  
  # The numbers in the final equation.
  A = 1
  B_num = 3
  B_den = 8
  C = 8
  D_num = 3
  D_den = 4
  
  # Format the equation string to explicitly show all the numbers.
  equation = (
      "The analytical solution for the normalized density profile rho'(xi) is:\n\n"
      "                  exp(({B_num}/{B_den}) * xi)\n"
      "rho'(xi) = {A} + --------------------------------------\n"
      "            sqrt({C} + exp(({D_num}/{D_den}) * xi))"
  ).format(
      A=A,
      B_num=B_num, B_den=B_den,
      C=C,
      D_num=D_num, D_den=D_den
  )

  # Let's print a more machine-readable version as well.
  equation_oneline = "rho'(xi) = {A} + exp(({B_num}/{B_den}) * xi) / sqrt({C} + exp(({D_num}/{D_den}) * xi))".format(
      A=A,
      B_num=B_num, B_den=B_den,
      C=C,
      D_num=D_num, D_den=D_den
  )

  print(equation_oneline)

if __name__ == "__main__":
  print_analytical_solution()