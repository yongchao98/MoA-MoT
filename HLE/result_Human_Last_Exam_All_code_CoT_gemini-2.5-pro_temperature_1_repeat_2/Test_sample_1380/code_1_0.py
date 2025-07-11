# The user is asking for the name of a physical process that Louis Bachelier
# found to be mathematically identical to the behavior of stock prices.
# This process is the random, erratic movement of particles suspended in a fluid.
# The script will store this name in a variable and then print it.

def get_bachelier_discovery():
  """
  This function identifies and returns the name of the physical process
  central to Louis Bachelier's work on financial markets.
  """
  process_name = "Brownian motion"
  return process_name

# Print the name of the discovered process
print(f"The physical process is: {get_bachelier_discovery()}")