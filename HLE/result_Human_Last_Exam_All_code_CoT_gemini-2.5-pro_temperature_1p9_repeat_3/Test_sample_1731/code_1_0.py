import sys

def solve_statistical_mechanics():
    """
    Derives and prints the equilibrium values of entropy and mean energy for a
    photon gas using principles of statistical mechanics and large deviation theory.
    """
    
    # Use a raw string for the explanation to handle backslashes and special chars.
    explanation = r"""
Derivation of Equilibrium Entropy and Mean Energy for a Photon Gas
=====================================================================
This script derives the equilibrium properties of a gas of light quanta (photons)
using principles from statistical mechanics. These principles, such as maximizing
entropy, are deeply connected to large deviation theorems.

*   The Boltzmann-Sanov theorem states that the probability of observing an empirical
    distribution is exponentially related to its "distance" (KL divergence) from the
    true distribution. Maximizing entropy is equivalent to finding the most probable
    macrostate, which minimizes this divergence.

*   The Cramer-Chernoff theorem relates probabilities of large deviations of a sum
    of random variables to a rate function. In physics, entropy acts as this rate
    function for energy, and it's related to the partition function via a
    Legendre transform.

Our plan follows these ideas:
1. Define the entropy (S) for a system of bosons (photons).
2. Maximize S with a fixed mean energy (E) to find the most probable, i.e.,
   equilibrium, distribution of photons across energy levels.
3. Use this equilibrium distribution to state the final equations for mean energy
   and entropy.

---------------------------------
Step 1: Entropy of a Bose Gas
---------------------------------
Let's define our terms:
 - epsilon_i : The energy of the i-th energy level.
 - g_i       : The number of states (degeneracy) at energy level epsilon_i.
 - N_i       : The number of photons occupying the energy level epsilon_i.
 - k_B       : The Boltzmann constant.
 - beta      : The inverse temperature, beta = 1 / (k_B * T).

The statistical entropy (S) for bosons, after applying Stirling's approximation, is:
S / k_B = sum_i [ (N_i + g_i)*ln(N_i + g_i) - N_i*ln(N_i) - g_i*ln(g_i) ]

-------------------------------------------
Step 2: Maximizing Entropy (Finding Equilibrium)
-------------------------------------------
We maximize S under the constraint of a fixed total energy E = sum_i (N_i * epsilon_i).
Since photons can be created and destroyed, their total number is not conserved (which
is equivalent to setting their chemical potential to zero).

Using the method of Lagrange multipliers, we maximize L = S/k_B - beta * (sum_i(N_i*epsilon_i) - E).
Setting the derivative dL/dN_i = 0 gives:
d(S/k_B)/dN_i = beta * epsilon_i
ln( (N_i + g_i) / N_i ) = beta * epsilon_i

Solving for N_i gives the equilibrium occupation number, denoted as <N_i>:
<N_i> = g_i / (exp(beta * epsilon_i) - 1)

This is the Bose-Einstein distribution for photons.

------------------------------------------
Step 3: Equilibrium Mean Energy and Entropy
------------------------------------------
Now we substitute the equilibrium distribution <N_i> back into the formulas for
the total energy and entropy to get their equilibrium values.

*** EQUILIBRIUM MEAN ENERGY (U) ***
The mean energy U is the sum of energies of all photons at equilibrium.

Final Equation:
U = sum_over_i [ (g_i * epsilon_i) / (exp(beta * epsilon_i) - 1) ]

where:
  U           = Equilibrium Mean Energy
  sum_over_i  = Sum over all energy levels i
  g_i         = Degeneracy of energy level i
  epsilon_i   = Energy of level i
  beta        = 1 / (k_B * T)
  exp()       = The exponential function

*** EQUILIBRIUM ENTROPY (S) ***
The entropy is found by substituting <N_i> back into the entropy formula or by using
the thermodynamic relation S = (U + P*V)/T = U/T + k_B * ln(Z), where Z is the partition
function. Both methods yield the same result.

Final Equation:
S / k_B = sum_over_i [ (g_i * beta * epsilon_i) / (exp(beta * epsilon_i) - 1) - g_i * ln(1 - exp(-beta * epsilon_i)) ]

where:
  S           = Equilibrium Entropy
  k_B         = Boltzmann constant
  ln()        = The natural logarithm
  All other symbols are as defined above.
"""
    print(explanation)

solve_statistical_mechanics()

# Final Answer Check. The question asks for the values. Since this is a theoretical derivation,
# the "values" are the final equilibrium equations. The code above prints them clearly.
# The user wants an answer in the format <<<answer content>>>. 
# It is unclear what should go here, as the result is two complex equations, not a single value.
# I will output the final derived formulas as the answer.

final_answer = "Mean Energy U = sum_i[g_i*epsilon_i / (exp(beta*epsilon_i) - 1)]; Entropy S/k_B = sum_i[g_i*beta*epsilon_i / (exp(beta*epsilon_i) - 1) - g_i*ln(1 - exp(-beta*epsilon_i))]"
sys.stdout.write(f"\n<<<{final_answer}>>>\n")
